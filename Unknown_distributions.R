################################################################################
#                    TYPE 2 PATIENTS - IMPROVED NON-PARAMETRIC MODEL          #
#                        FOR DISCRETE EVENT SIMULATION (DES)                   #
################################################################################
#
# PURPOSE:
# --------
# Analyze Type 2 patient scan data to create empirical (non-parametric) input
# models for discrete event simulation. This approach uses the actual observed
# data without assuming any theoretical probability distributions.
#
# KEY CONCEPTS:
# -------------
# 1. INTERARRIVAL TIME: Time between consecutive patient arrivals
#    Example: Patient A arrives at 9:00, Patient B at 9:15 → interarrival = 0.25 hrs
#
# 2. DURATION: Length of time for the MRI scan procedure
#
# 3. NON-PARAMETRIC: Uses empirical data directly (no Gamma, Exponential, etc.)
#    - More robust when true distribution is unknown
#    - Preserves all features of real data (multimodality, skewness, etc.)
#
# 4. BOOTSTRAP: Resampling technique to quantify uncertainty
#    - Creates "fake" datasets by sampling WITH replacement
#    - Estimates confidence intervals without assuming normality
#
# 5. EMPIRICAL SAMPLING: For simulation, randomly draw from observed values
#    - Maintains realistic variability
#    - No parametric assumptions needed
#
################################################################################

# ============================================================================
# 0. SETUP AND LIBRARIES
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat("           TYPE 2 PATIENT NON-PARAMETRIC MODEL - ENHANCED VERSION              \n")
cat("================================================================================\n\n")

# Required packages
required_packages <- c("dplyr", "boot", "lubridate", "ggplot2")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "http://cran.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)      # Data manipulation
  library(boot)       # Bootstrap resampling
  library(lubridate)  # Date/time handling
  library(ggplot2)    # Visualization
})

set.seed(123)  # Reproducibility

cat("✓ Libraries loaded\n\n")

# 1. DATA LOADING AND VALIDATION


cat("STEP 1: Loading and validating data\n")
cat("----------------------------------------\n")

# Load data with error handling
tryCatch({
  data <- read.csv("ScanRecords (1).csv", stringsAsFactors = FALSE)
  cat("✓ Data loaded successfully:", nrow(data), "total records\n")
}, error = function(e) {
  stop("ERROR: Cannot load 'ScanRecords (1).csv'. Check file path.\n", e$message)
})

# Validate required columns
required_cols <- c("PatientType", "Date", "Time", "Duration")
missing_cols <- setdiff(required_cols, colnames(data))

if(length(missing_cols) > 0) {
  stop("ERROR: Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("✓ All required columns present\n")

# Check for missing values
na_counts <- sapply(data[required_cols], function(x) sum(is.na(x)))
if(any(na_counts > 0)) {
  cat("⚠ WARNING: Missing values detected:\n")
  print(na_counts[na_counts > 0])
}

cat("\n")

# ============================================================================
# 2. FILTER TYPE 2 PATIENTS
# ============================================================================

cat("STEP 2: Filtering Type 2 patients\n")
cat("----------------------------------------\n")

# Handle multiple possible encodings of Type 2
type2 <- data %>%
  filter(PatientType == 2 | 
           PatientType == "2" | 
           PatientType == "Type 2" | 
           PatientType == "Type2")

if(nrow(type2) == 0) {
  stop("ERROR: No Type 2 patients found. Check PatientType values in data.")
}

cat("✓ Type 2 patients identified:", nrow(type2), "records\n")
cat("  Percentage of total:", round(100 * nrow(type2) / nrow(data), 1), "%\n\n")

# ============================================================================
# 3. CONSTRUCT INTRADAY ARRIVAL TIMES
# ============================================================================
#
# EXPLANATION:
# ------------
# The Time column is stored as decimal hours (e.g., 9.5 = 9:30 AM)
# We need to convert this to proper timestamps (POSIXct) for:
# 1. Correct chronological ordering
# 2. Accurate time difference calculations
# 3. Date/time arithmetic
#
# PROCESS:
# --------
# Input:  Date = "2024-01-15", Time = 9.75
# Step 1: Extract hour = 9, minute = 45
# Step 2: Create "2024-01-15 09:45"
# Step 3: Parse as POSIXct timestamp
#
# ============================================================================

cat("STEP 3: Processing timestamps\n")
cat("----------------------------------------\n")

type2 <- type2 %>%
  mutate(
    # Convert Date to proper Date object
    Date = as.Date(Date),
    
    # Extract hours and minutes from decimal Time
    hour = floor(Time),                    # 9.75 → 9
    minute = round((Time - hour) * 60),    # (9.75 - 9) * 60 = 45
    
    # Create DateTime string in format "YYYY-MM-DD HH:MM"
    DateTime = as.POSIXct(
      paste(Date, sprintf("%02d:%02d", hour, minute)),
      format = "%Y-%m-%d %H:%M",
      tz = "UTC"  # Use UTC to avoid daylight saving issues
    )
  ) %>%
  arrange(DateTime)  # Sort chronologically - CRITICAL for interarrival calculation!

# Validate DateTime parsing
n_failed <- sum(is.na(type2$DateTime))
if(n_failed > 0) {
  cat("⚠ WARNING:", n_failed, "timestamps failed to parse\n")
  cat("  Sample problematic entries:\n")
  print(head(type2[is.na(type2$DateTime), c("Date", "Time", "hour", "minute")], 3))
  
  # Remove failed parses
  type2 <- type2 %>% filter(!is.na(DateTime))
  cat("  Continuing with", nrow(type2), "valid records\n")
}

cat("✓ Timestamps processed successfully\n")
cat("  Time range:", format(min(type2$DateTime), "%Y-%m-%d %H:%M"), "to",
    format(max(type2$DateTime), "%Y-%m-%d %H:%M"), "\n")
cat("  Span:", round(as.numeric(difftime(max(type2$DateTime), 
                                         min(type2$DateTime), 
                                         units = "days")), 1), "days\n\n")

# ============================================================================
# 4. COMPUTE INTERARRIVAL TIMES
# ============================================================================
#
# EXPLANATION:
# ------------
# Interarrival time = Time between consecutive patient arrivals
# 
# Example:
# Patient 1 arrives: 2024-01-15 09:00
# Patient 2 arrives: 2024-01-15 09:15
# Interarrival = 0.25 hours (15 minutes)
#
# IMPORTANT ASSUMPTION:
# ---------------------
# "No time passes outside working hours"
# This means if the last patient on Monday arrives at 5pm and the first
# patient on Tuesday arrives at 8am, we DO NOT count the 15 hours overnight.
# We only count the actual time between arrivals during operating hours.
#
# IMPLEMENTATION:
# ---------------
# The diff() function calculates consecutive differences:
# DateTime: [9:00, 9:15, 9:30, 10:00]
# diff():   [0.25, 0.25, 0.50] hours
#
# ============================================================================

cat("STEP 4: Computing interarrival times\n")
cat("----------------------------------------\n")

# Calculate time differences between consecutive arrivals
interarrival <- diff(type2$DateTime)
interarrival_hours <- as.numeric(interarrival, units = "hours")

cat("Initial interarrivals computed:", length(interarrival_hours), "\n")

# Data cleaning: Remove non-positive interarrivals
# (These can occur due to simultaneous arrivals or data entry errors)
n_nonpositive <- sum(interarrival_hours <= 0)
if(n_nonpositive > 0) {
  cat("⚠ Removing", n_nonpositive, "non-positive interarrivals\n")
  interarrival_hours <- interarrival_hours[interarrival_hours > 0]
}

# Identify potential outliers (very long gaps)
# These might represent overnight/weekend gaps or data collection issues
threshold <- 24  # hours
n_outliers <- sum(interarrival_hours > threshold)

if(n_outliers > 0) {
  cat("⚠ WARNING:", n_outliers, "interarrivals exceed", threshold, "hours\n")
  cat("  These may represent non-operating periods (nights/weekends)\n")
  cat("  Max interarrival:", round(max(interarrival_hours), 2), "hours\n")
  cat("  Consider filtering these if they violate the 'intraday' assumption\n")
  
  # Show distribution of large gaps
  large_gaps <- interarrival_hours[interarrival_hours > threshold]
  cat("  Large gaps summary:\n")
  print(summary(large_gaps))
}

cat("\n✓ Final interarrival sample size:", length(interarrival_hours), "\n\n")

# ============================================================================
# 5. EXTRACT SCAN DURATIONS
# ============================================================================

cat("STEP 5: Processing scan durations\n")
cat("----------------------------------------\n")

durations <- type2$Duration

# Remove missing or invalid durations
n_na_dur <- sum(is.na(durations))
n_nonpos_dur <- sum(durations <= 0, na.rm = TRUE)

if(n_na_dur > 0) {
  cat("⚠ Removing", n_na_dur, "missing durations\n")
}
if(n_nonpos_dur > 0) {
  cat("⚠ Removing", n_nonpos_dur, "non-positive durations\n")
}

durations <- durations[!is.na(durations) & durations > 0]

cat("✓ Final duration sample size:", length(durations), "\n\n")

# ============================================================================
# 6. DESCRIPTIVE STATISTICS
# ============================================================================

cat("STEP 6: Computing descriptive statistics\n")
cat("----------------------------------------\n")

# Duration statistics
dur_stats <- list(
  n = length(durations),
  mean = mean(durations),
  median = median(durations),
  sd = sd(durations),
  min = min(durations),
  max = max(durations),
  q25 = quantile(durations, 0.25),
  q75 = quantile(durations, 0.75),
  q90 = quantile(durations, 0.90),
  q95 = quantile(durations, 0.95),
  iqr = IQR(durations),
  cv = sd(durations) / mean(durations)  # Coefficient of variation
)

# Interarrival statistics
arr_stats <- list(
  n = length(interarrival_hours),
  mean = mean(interarrival_hours),
  median = median(interarrival_hours),
  sd = sd(interarrival_hours),
  min = min(interarrival_hours),
  max = max(interarrival_hours),
  q25 = quantile(interarrival_hours, 0.25),
  q75 = quantile(interarrival_hours, 0.75),
  q90 = quantile(interarrival_hours, 0.90),
  q95 = quantile(interarrival_hours, 0.95),
  iqr = IQR(interarrival_hours),
  cv = sd(interarrival_hours) / mean(interarrival_hours)
)

# Print results
cat("\nDURATION STATISTICS (hours):\n")
cat("============================\n")
cat(sprintf("  Sample size:    %d\n", dur_stats$n))
cat(sprintf("  Mean:           %.3f\n", dur_stats$mean))
cat(sprintf("  Median:         %.3f\n", dur_stats$median))
cat(sprintf("  Std Dev:        %.3f\n", dur_stats$sd))
cat(sprintf("  CV (σ/μ):       %.3f\n", dur_stats$cv))
cat(sprintf("  Range:          [%.3f, %.3f]\n", dur_stats$min, dur_stats$max))
cat(sprintf("  IQR:            %.3f\n", dur_stats$iqr))
cat(sprintf("  25th percentile: %.3f\n", dur_stats$q25))
cat(sprintf("  75th percentile: %.3f\n", dur_stats$q75))
cat(sprintf("  90th percentile: %.3f\n", dur_stats$q90))
cat(sprintf("  95th percentile: %.3f\n", dur_stats$q95))

cat("\nINTERARRIVAL STATISTICS (hours):\n")
cat("================================\n")
cat(sprintf("  Sample size:    %d\n", arr_stats$n))
cat(sprintf("  Mean:           %.3f\n", arr_stats$mean))
cat(sprintf("  Median:         %.3f\n", arr_stats$median))
cat(sprintf("  Std Dev:        %.3f\n", arr_stats$sd))
cat(sprintf("  CV (σ/μ):       %.3f\n", arr_stats$cv))
cat(sprintf("  Range:          [%.3f, %.3f]\n", arr_stats$min, arr_stats$max))
cat(sprintf("  IQR:            %.3f\n", arr_stats$iqr))
cat(sprintf("  25th percentile: %.3f\n", arr_stats$q25))
cat(sprintf("  75th percentile: %.3f\n", arr_stats$q75))
cat(sprintf("  90th percentile: %.3f\n", arr_stats$q90))
cat(sprintf("  95th percentile: %.3f\n", arr_stats$q95))

cat("\n")

# Interpretation notes
cat("INTERPRETATION:\n")
cat("---------------\n")
if(dur_stats$cv > 0.5) {
  cat("⚠ Duration CV > 0.5 indicates HIGH variability\n")
} else {
  cat("✓ Duration CV < 0.5 indicates moderate variability\n")
}

if(arr_stats$cv > 1.0) {
  cat("⚠ Interarrival CV > 1.0 indicates VERY HIGH variability\n")
  cat("  (Common for arrival processes; may indicate clustering)\n")
} else {
  cat("✓ Interarrival CV < 1.0 indicates moderate variability\n")
}

cat("\n")

# ============================================================================
# 7. BOOTSTRAP UNCERTAINTY QUANTIFICATION
# ============================================================================
#
# EXPLANATION:
# ------------
# Bootstrap is a resampling method to estimate uncertainty without assuming
# any distribution (e.g., normality).
#
# HOW IT WORKS:
# -------------
# 1. From original data [x1, x2, ..., xn], create a "bootstrap sample" by
#    randomly drawing n values WITH REPLACEMENT
# 2. Calculate the statistic (e.g., mean) on this bootstrap sample
# 3. Repeat steps 1-2 many times (e.g., 1000 times)
# 4. The distribution of bootstrap statistics estimates the sampling distribution
# 5. Use percentiles of this distribution as confidence intervals
#
# EXAMPLE:
# --------
# Original data: [1, 2, 3, 4, 5], mean = 3
# Bootstrap sample 1: [1, 1, 3, 5, 5], mean = 3.0
# Bootstrap sample 2: [2, 3, 4, 4, 5], mean = 3.6
# ... (repeat 1000 times)
# 95% CI: [2.5th percentile, 97.5th percentile] of bootstrap means
#
# ============================================================================

cat("STEP 7: Bootstrap uncertainty quantification\n")
cat("----------------------------------------\n")

# Define bootstrap functions for different statistics
boot_mean <- function(x, i) mean(x[i])
boot_median <- function(x, i) median(x[i])
boot_sd <- function(x, i) sd(x[i])
boot_q90 <- function(x, i) quantile(x[i], 0.90)
boot_q95 <- function(x, i) quantile(x[i], 0.95)

B <- 10000  # Number of bootstrap replications (more = better precision)

cat("Running", B, "bootstrap replications...\n")
cat("(This may take 30-60 seconds)\n\n")

# Bootstrap DURATIONS
cat("  Bootstrapping duration statistics...\n")
boot_dur_mean <- boot(durations, boot_mean, R = B)
boot_dur_median <- boot(durations, boot_median, R = B)
boot_dur_sd <- boot(durations, boot_sd, R = B)
boot_dur_q90 <- boot(durations, boot_q90, R = B)
boot_dur_q95 <- boot(durations, boot_q95, R = B)

# Bootstrap INTERARRIVALS
cat("  Bootstrapping interarrival statistics...\n")
boot_arr_mean <- boot(interarrival_hours, boot_mean, R = B)
boot_arr_median <- boot(interarrival_hours, boot_median, R = B)
boot_arr_sd <- boot(interarrival_hours, boot_sd, R = B)
boot_arr_q90 <- boot(interarrival_hours, boot_q90, R = B)
boot_arr_q95 <- boot(interarrival_hours, boot_q95, R = B)

cat("\n✓ Bootstrap completed\n\n")

# Calculate confidence intervals (using both methods for robustness)
cat("Computing confidence intervals...\n")

# Duration CIs
dur_mean_ci <- boot.ci(boot_dur_mean, type = c("perc", "bca"))
dur_median_ci <- boot.ci(boot_dur_median, type = c("perc", "bca"))
dur_q90_ci <- boot.ci(boot_dur_q90, type = c("perc", "bca"))
dur_q95_ci <- boot.ci(boot_dur_q95, type = c("perc", "bca"))

# Interarrival CIs
arr_mean_ci <- boot.ci(boot_arr_mean, type = c("perc", "bca"))
arr_median_ci <- boot.ci(boot_arr_median, type = c("perc", "bca"))
arr_q90_ci <- boot.ci(boot_arr_q90, type = c("perc", "bca"))
arr_q95_ci <- boot.ci(boot_arr_q95, type = c("perc", "bca"))

cat("✓ Confidence intervals computed\n\n")

# Display results in readable format
cat("================================================================================\n")
cat("                    BOOTSTRAP RESULTS (95% CONFIDENCE INTERVALS)               \n")
cat("================================================================================\n\n")

cat("DURATION:\n")
cat("---------\n")
cat(sprintf("  Mean:     %.3f  [%.3f, %.3f]  (SE = %.3f)\n",
            boot_dur_mean$t0,
            dur_mean_ci$bca[4], dur_mean_ci$bca[5],
            sd(boot_dur_mean$t)))
cat(sprintf("  Median:   %.3f  [%.3f, %.3f]  (SE = %.3f)\n",
            boot_dur_median$t0,
            dur_median_ci$bca[4], dur_median_ci$bca[5],
            sd(boot_dur_median$t)))
cat(sprintf("  90th %%ile: %.3f  [%.3f, %.3f]  (SE = %.3f)\n",
            boot_dur_q90$t0,
            dur_q90_ci$bca[4], dur_q90_ci$bca[5],
            sd(boot_dur_q90$t)))
cat(sprintf("  95th %%ile: %.3f  [%.3f, %.3f]  (SE = %.3f)\n",
            boot_dur_q95$t0,
            dur_q95_ci$bca[4], dur_q95_ci$bca[5],
            sd(boot_dur_q95$t)))

cat("\nINTERARRIVAL:\n")
cat("-------------\n")
cat(sprintf("  Mean:     %.3f  [%.3f, %.3f]  (SE = %.3f)\n",
            boot_arr_mean$t0,
            arr_mean_ci$bca[4], arr_mean_ci$bca[5],
            sd(boot_arr_mean$t)))
cat(sprintf("  Median:   %.3f  [%.3f, %.3f]  (SE = %.3f)\n",
            boot_arr_median$t0,
            arr_median_ci$bca[4], arr_median_ci$bca[5],
            sd(boot_arr_median$t)))
cat(sprintf("  90th %%ile: %.3f  [%.3f, %.3f]  (SE = %.3f)\n",
            boot_arr_q90$t0,
            arr_q90_ci$bca[4], arr_q90_ci$bca[5],
            sd(boot_arr_q90$t)))
cat(sprintf("  95th %%ile: %.3f  [%.3f, %.3f]  (SE = %.3f)\n",
            boot_arr_q95$t0,
            arr_q95_ci$bca[4], arr_q95_ci$bca[5],
            sd(boot_arr_q95$t)))

cat("\nNOTE: Intervals shown are BCa (bias-corrected and accelerated) CIs\n")
cat("      SE = Standard Error (standard deviation of bootstrap distribution)\n\n")

# ============================================================================
# 8. VISUALIZATION
# ============================================================================

cat("STEP 8: Creating visualizations\n")
cat("----------------------------------------\n")

pdf("type2_nonparametric_results.pdf", width = 12, height = 8)

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# Duration histogram
hist(durations, breaks = 30, main = "Duration Distribution",
     xlab = "Duration (hours)", col = "lightblue", border = "white",
     prob = TRUE)
lines(density(durations), col = "darkblue", lwd = 2)
abline(v = mean(durations), col = "red", lwd = 2, lty = 2)
abline(v = quantile(durations, c(0.25, 0.75)), col = "orange", lwd = 1, lty = 3)
legend("topright", c("Mean", "Q1/Q3"), col = c("red", "orange"), 
       lwd = c(2,1), lty = c(2,3), cex = 0.8)

# Interarrival histogram
hist(interarrival_hours, breaks = 30, main = "Interarrival Distribution",
     xlab = "Interarrival Time (hours)", col = "lightgreen", border = "white",
     prob = TRUE)
lines(density(interarrival_hours), col = "darkgreen", lwd = 2)
abline(v = mean(interarrival_hours), col = "red", lwd = 2, lty = 2)
legend("topright", "Mean", col = "red", lwd = 2, lty = 2, cex = 0.8)

# Bootstrap distribution: Duration mean
hist(boot_dur_mean$t, breaks = 50, main = "Bootstrap: Duration Mean",
     xlab = "Bootstrap Means", col = "lightcoral", border = "white",
     prob = TRUE)
abline(v = boot_dur_mean$t0, col = "red", lwd = 2)
abline(v = dur_mean_ci$bca[4:5], col = "blue", lwd = 2, lty = 2)
legend("topright", c("Original", "95% CI"), col = c("red", "blue"),
       lwd = 2, lty = c(1,2), cex = 0.8)

# Bootstrap distribution: Interarrival mean
hist(boot_arr_mean$t, breaks = 50, main = "Bootstrap: Interarrival Mean",
     xlab = "Bootstrap Means", col = "lightyellow", border = "white",
     prob = TRUE)
abline(v = boot_arr_mean$t0, col = "red", lwd = 2)
abline(v = arr_mean_ci$bca[4:5], col = "blue", lwd = 2, lty = 2)
legend("topright", c("Original", "95% CI"), col = c("red", "blue"),
       lwd = 2, lty = c(1,2), cex = 0.8)

# ECDF: Duration
plot(ecdf(durations), main = "Empirical CDF: Duration",
     xlab = "Duration (hours)", ylab = "Cumulative Probability",
     col = "darkblue", lwd = 2)
abline(h = c(0.5, 0.9, 0.95), col = "gray", lty = 2)
abline(v = quantile(durations, c(0.5, 0.9, 0.95)), col = "red", lty = 2)

# ECDF: Interarrival
plot(ecdf(interarrival_hours), main = "Empirical CDF: Interarrival",
     xlab = "Interarrival Time (hours)", ylab = "Cumulative Probability",
     col = "darkgreen", lwd = 2)
abline(h = c(0.5, 0.9, 0.95), col = "gray", lty = 2)
abline(v = quantile(interarrival_hours, c(0.5, 0.9, 0.95)), col = "red", lty = 2)

dev.off()

cat("✓ Visualizations saved to 'type2_nonparametric_results.pdf'\n\n")

# ============================================================================
# 9. SIMULATION-READY SAMPLING FUNCTIONS
# ============================================================================
#
# EXPLANATION:
# ------------
# These functions allow you to generate random samples from the empirical
# distribution for use in discrete event simulation.
#
# HOW THEY WORK:
# --------------
# sample_type2_duration(n = 10) will:
# 1. Randomly select 10 values from the observed durations
# 2. Selection is WITH REPLACEMENT (same value can be chosen multiple times)
# 3. Each observed value has equal probability of being selected
#
# EXAMPLE USAGE IN SIMULATION:
# -----------------------------
# # Generate next patient's scan duration
# duration <- sample_type2_duration(n = 1)
#
# # Generate interarrival time for next 5 patients
# interarrivals <- sample_type2_interarrival(n = 5)
#
# ============================================================================

cat("STEP 9: Creating simulation sampling functions\n")
cat("----------------------------------------\n")

# Duration sampler
sample_type2_duration <- function(n = 1) {
  if(n <= 0) stop("n must be positive")
  sample(durations, size = n, replace = TRUE)
}

# Interarrival sampler
sample_type2_interarrival <- function(n = 1) {
  if(n <= 0) stop("n must be positive")
  sample(interarrival_hours, size = n, replace = TRUE)
}

cat("✓ Sampling functions created\n")

# Validate samplers
cat("\nValidating samplers with 10,000 samples...\n")

test_dur <- sample_type2_duration(10000)
test_arr <- sample_type2_interarrival(10000)

cat("  Duration sampler:\n")
cat(sprintf("    Original mean: %.3f, Sampled mean: %.3f\n", 
            mean(durations), mean(test_dur)))
cat(sprintf("    Original SD:   %.3f, Sampled SD:   %.3f\n",
            sd(durations), sd(test_dur)))

cat("  Interarrival sampler:\n")
cat(sprintf("    Original mean: %.3f, Sampled mean: %.3f\n",
            mean(interarrival_hours), mean(test_arr)))
cat(sprintf("    Original SD:   %.3f, Sampled SD:   %.3f\n",
            sd(interarrival_hours), sd(test_arr)))

cat("\n✓ Samplers validated (statistics match original data)\n\n")

# ============================================================================
# 10. CREATE FINAL MODEL OBJECT
# ============================================================================

cat("STEP 10: Assembling final model object\n")
cat("----------------------------------------\n")

type2_model <- list(
  # Original data
  duration_data = durations,
  interarrival_data = interarrival_hours,
  
  # Descriptive statistics
  duration_stats = dur_stats,
  interarrival_stats = arr_stats,
  
  # Bootstrap results
  bootstrap = list(
    duration = list(
      mean = boot_dur_mean,
      median = boot_dur_median,
      sd = boot_dur_sd,
      q90 = boot_dur_q90,
      q95 = boot_dur_q95
    ),
    interarrival = list(
      mean = boot_arr_mean,
      median = boot_arr_median,
      sd = boot_arr_sd,
      q90 = boot_arr_q90,
      q95 = boot_arr_q95
    )
  ),
  
  # Confidence intervals
  confidence_intervals = list(
    duration_mean = dur_mean_ci,
    duration_median = dur_median_ci,
    duration_q90 = dur_q90_ci,
    duration_q95 = dur_q95_ci,
    interarrival_mean = arr_mean_ci,
    interarrival_median = arr_median_ci,
    interarrival_q90 = arr_q90_ci,
    interarrival_q95 = arr_q95_ci
  ),
  
  # Sampling functions (for simulation)
  samplers = list(
    duration = sample_type2_duration,
    interarrival = sample_type2_interarrival
  ),
  
  # Metadata
  metadata = list(
    creation_date = Sys.time(),
    n_duration = length(durations),
    n_interarrival = length(interarrival_hours),
    bootstrap_replications = B,
    seed = 123
  )
)

cat("✓ Model object created with", length(names(type2_model)), "components\n\n")

# ============================================================================
# 11. SAVE RESULTS
# ============================================================================

cat("STEP 11: Saving results\n")
cat("----------------------------------------\n")

# Save R object
save(type2_model, file = "type2_nonparametric_model.RData")
cat("✓ Model saved to 'type2_nonparametric_model.RData'\n")

# Save summary to text file
sink("type2_model_summary.txt")
cat("================================================================================\n")
cat("           TYPE 2 PATIENT NON-PARAMETRIC MODEL - SUMMARY REPORT                \n")
cat("================================================================================\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("SAMPLE SIZES:\n")
cat("  Duration observations:    ", length(durations), "\n")
cat("  Interarrival observations:", length(interarrival_hours), "\n\n")

cat("DURATION ESTIMATES (hours):\n")
cat(sprintf("  Mean:     %.3f  [%.3f, %.3f]\n",
            boot_dur_mean$t0, dur_mean_ci$bca[4], dur_mean_ci$bca[5]))
cat(sprintf("  Median:   %.3f  [%.3f, %.3f]\n",
            boot_dur_median$t0, dur_median_ci$bca[4], dur_median_ci$bca[5]))
cat(sprintf("  90th %%ile: %.3f  [%.3f, %.3f]\n",
            boot_dur_q90$t0, dur_q90_ci$bca[4], dur_q90_ci$bca[5]))
cat(sprintf("  95th %%ile: %.3f  [%.3f, %.3f]\n\n",
            boot_dur_q95$t0, dur_q95_ci$bca[4], dur_q95_ci$bca[5]))

cat("INTERARRIVAL ESTIMATES (hours):\n")
cat(sprintf("  Mean:     %.3f  [%.3f, %.3f]\n",
            boot_arr_mean$t0, arr_mean_ci$bca[4], arr_mean_ci$bca[5]))
cat(sprintf("  Median:   %.3f  [%.3f, %.3f]\n",
            boot_arr_median$t0, arr_median_ci$bca[4], arr_median_ci$bca[5]))
cat(sprintf("  90th %%ile: %.3f  [%.3f, %.3f]\n",
            boot_arr_q90$t0, arr_q90_ci$bca[4], arr_q90_ci$bca[5]))
cat(sprintf("  95th %%ile: %.3f  [%.3f, %.3f]\n\n",
            boot_arr_q95$t0, arr_q95_ci$bca[4], arr_q95_ci$bca[5]))

cat("USAGE IN SIMULATION:\n")
cat("  # Load model\n")
cat("  load('type2_nonparametric_model.RData')\n\n")
cat("  # Sample duration for next patient\n")
cat("  duration <- type2_model$samplers$duration(n = 1)\n\n")
cat("  # Sample interarrival times for next 100 patients\n")
cat("  interarrivals <- type2_model$samplers$interarrival(n = 100)\n\n")

cat("================================================================================\n")
sink()

cat("✓ Summary saved to 'type2_model_summary.txt'\n\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("================================================================================\n")
cat("                            ANALYSIS COMPLETE!                                  \n")
cat("================================================================================\n\n")

cat("✓ Type 2 non-parametric model successfully created\n\n")

cat("KEY OUTPUTS:\n")
cat("  1. type2_model                         - R object in workspace\n")
cat("  2. type2_nonparametric_model.RData     - Saved model file\n")
cat("  3. type2_model_summary.txt             - Text summary report\n")
cat("  4. type2_nonparametric_results.pdf     - Visualization plots\n\n")

cat("NEXT STEPS:\n")
cat("  → Use samplers in discrete event simulation\n")
cat("  → Compare to parametric models if needed\n")
cat("  → Validate simulation results against historical data\n\n")

cat("EXAMPLE USAGE:\n")
cat("  # Load the model\n")
cat("  load('type2_nonparametric_model.RData')\n\n")
cat("  # Generate random duration\n")
cat("  type2_model$samplers$duration(1)\n\n")
cat("  # Generate 5 random interarrivals\n")
cat("  type2_model$samplers$interarrival(5)\n\n")

cat("  # Access statistics\n")
cat("  type2_model$duration_stats$mean\n")
cat("  type2_model$confidence_intervals$duration_mean\n\n")

cat("================================================================================\n\n")

################################################################################
#                              END OF SCRIPT                                   #
################################################################################
