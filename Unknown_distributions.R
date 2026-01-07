library(readr)
library(ggplot2)

# Read data
data <- read_csv("ScanRecords (1).csv")

# Count patients correctly
numberdays <- length(unique(data$Date))
numberPatients <- matrix(data = 0, nrow=numberdays, ncol = 2)
colnames(numberPatients) <- c("Type1","Type2")
counter <- 1

for(i in 1:length(data$Date)){
  if(i > 1 && data$Date[i] != data$Date[i-1]){
    counter = counter + 1
  }
  if(data$PatientType[i] == "Type 2"){
    numberPatients[counter,2] <- numberPatients[counter,2] + 1  # Fixed: Type 2 goes to column 2
  }
  else{
    numberPatients[counter,1] <- numberPatients[counter,1] + 1  # Type 1 goes to column 1
  }
}

# Extract Type 2 counts
type2_counts <- numberPatients[,2]

# Summary statistics
cat("Summary Statistics for Type 2 patients per day:\n")
print(summary(type2_counts))
cat("\nStandard Deviation:", sd(type2_counts), "\n")
cat("Variance:", var(type2_counts), "\n")

# Visualizations
par(mfrow=c(2,3))

# 1. Histogram with density overlay
hist(type2_counts, breaks=15, probability=TRUE, 
     main="Type 2 Patients per Day\nHistogram with Density", 
     xlab="Number of Type 2 Patients", col="lightblue", border="black")
lines(density(type2_counts), col="red", lwd=2)

# 2. Boxplot
boxplot(type2_counts, main="Type 2 Patients per Day\nBoxplot", 
        ylab="Number of Patients", col="lightgreen")

# 3. Q-Q plot for normal distribution
qqnorm(type2_counts, main="Q-Q Plot (Normal)")
qqline(type2_counts, col="red", lwd=2)

# 4. Bar plot of frequency
barplot(table(type2_counts), main="Frequency of Daily Counts", 
        xlab="Number of Type 2 Patients", ylab="Frequency", col="coral")

# 5. Time series plot
plot(type2_counts, type="l", main="Type 2 Patients Over Time", 
     xlab="Day", ylab="Number of Patients", col="blue")
points(type2_counts, pch=19, col="blue", cex=0.5)

# 6. ECDF (Empirical Cumulative Distribution Function)
plot(ecdf(type2_counts), main="Empirical CDF", 
     xlab="Number of Type 2 Patients", ylab="Cumulative Probability")

par(mfrow=c(1,1))

# Kolmogorov-Smirnov Tests
cat("\n=== KOLMOGOROV-SMIRNOV TESTS ===\n\n")

# Test 1: Normal Distribution
mean_type2 <- mean(type2_counts)
sd_type2 <- sd(type2_counts)
ks_normal <- ks.test(type2_counts, "pnorm", mean_type2, sd_type2)
cat("1. Normal Distribution Test:\n")
print(ks_normal)
cat("\n")

# Test 2: Poisson Distribution
lambda_type2 <- mean(type2_counts)
ks_poisson <- ks.test(type2_counts, "ppois", lambda_type2)
cat("2. Poisson Distribution Test:\n")
print(ks_poisson)
cat("\n")

# Test 3: Negative Binomial (if variance > mean, indicates overdispersion)
if(var(type2_counts) > mean(type2_counts)){
  cat("Note: Variance > Mean, suggesting possible Negative Binomial distribution\n")
  cat("Variance/Mean ratio:", var(type2_counts)/mean(type2_counts), "\n\n")
}

# Interpretation helper
cat("=== INTERPRETATION ===\n")
cat("If p-value > 0.05: Data is consistent with the tested distribution\n")
cat("If p-value < 0.05: Data significantly differs from the tested distribution\n\n")

# Better visualization comparing to theoretical distributions
par(mfrow=c(1,2))

# Compare to Normal
hist(type2_counts, breaks=15, probability=TRUE, 
     main="Type 2 vs Normal Distribution", 
     xlab="Number of Type 2 Patients", col="lightblue", border="black")
curve(dnorm(x, mean=mean_type2, sd=sd_type2), add=TRUE, col="red", lwd=2)
legend("topright", legend=c("Observed", "Normal"), 
       col=c("lightblue", "red"), lwd=c(10,2))

# Compare to Poisson
hist(type2_counts, breaks=15, probability=TRUE, 
     main="Type 2 vs Poisson Distribution", 
     xlab="Number of Type 2 Patients", col="lightgreen", border="black")
x_vals <- min(type2_counts):max(type2_counts)
points(x_vals, dpois(x_vals, lambda=lambda_type2), col="blue", pch=19, cex=1.5)
lines(x_vals, dpois(x_vals, lambda=lambda_type2), col="blue", lwd=2)
legend("topright", legend=c("Observed", "Poisson"), 
       col=c("lightgreen", "blue"), lwd=c(10,2))

par(mfrow=c(1,1))



library(readr)
library(ggplot2)
library(gridExtra)

# Read and process data
data <- read_csv("ScanRecords (1).csv")

# Count patients correctly
numberdays <- length(unique(data$Date))
numberPatients <- matrix(data = 0, nrow=numberdays, ncol = 2)
colnames(numberPatients) <- c("Type1","Type2")
counter <- 1

for(i in 1:length(data$Date)){
  if(i > 1 && data$Date[i] != data$Date[i-1]){
    counter = counter + 1
  }
  if(data$PatientType[i] == "Type 2"){
    numberPatients[counter,2] <- numberPatients[counter,2] + 1
  }
  else{
    numberPatients[counter,1] <- numberPatients[counter,1] + 1
  }
}

# Extract Type 2 counts
type2_counts <- numberPatients[,2]

cat("========================================\n")
cat("ORIGINAL DATA SUMMARY\n")
cat("========================================\n")
print(summary(type2_counts))
cat("Standard Deviation:", sd(type2_counts), "\n")
cat("Variance:", var(type2_counts), "\n")
cat("Sample Size:", length(type2_counts), "\n\n")

# ========================================
# PART 1: CHOOSING BOOTSTRAP METHOD
# ========================================

cat("========================================\n")
cat("BOOTSTRAP METHOD SELECTION\n")
cat("========================================\n\n")

# Assess distribution characteristics
mean_val <- mean(type2_counts)
var_val <- var(type2_counts)
skewness <- mean((type2_counts - mean_val)^3) / sd(type2_counts)^3

cat("Distribution Characteristics:\n")
cat("Mean:", mean_val, "\n")
cat("Variance:", var_val, "\n")
cat("Variance/Mean ratio:", var_val/mean_val, "\n")
cat("Skewness:", skewness, "\n\n")

# Test for normality and Poisson
shapiro_test <- shapiro.test(type2_counts)
ks_normal <- ks.test(type2_counts, "pnorm", mean_val, sd(type2_counts))
ks_poisson <- ks.test(type2_counts, "ppois", mean_val)

cat("Distribution Tests:\n")
cat("Shapiro-Wilk (normality): p-value =", shapiro_test$p.value, "\n")
cat("KS test (normal): p-value =", ks_normal$p.value, "\n")
cat("KS test (Poisson): p-value =", ks_poisson$p.value, "\n\n")

# Decision criteria
cat("BOOTSTRAP METHOD DECISION:\n")
cat("----------------------------------------\n")

use_parametric <- FALSE
justification <- ""

if(shapiro_test$p.value > 0.05 && ks_normal$p.value > 0.05){
  use_parametric <- TRUE
  justification <- "PARAMETRIC (Normal) Bootstrap chosen because:\n  - Data passes Shapiro-Wilk normality test (p > 0.05)\n  - Data fits normal distribution in KS test (p > 0.05)\n  - Parametric bootstrap is more efficient when distributional assumptions hold"
} else if(ks_poisson$p.value > 0.05){
  use_parametric <- TRUE
  justification <- "PARAMETRIC (Poisson) Bootstrap chosen because:\n  - Data fits Poisson distribution in KS test (p > 0.05)\n  - Count data naturally follows Poisson process\n  - Variance/Mean ratio close to 1 supports Poisson assumption"
} else {
  use_parametric <- FALSE
  justification <- "NON-PARAMETRIC Bootstrap chosen because:\n  - Data does not clearly fit standard parametric distributions\n  - Non-parametric method makes no distributional assumptions\n  - More robust to distributional misspecification\n  - Appropriate for unknown or complex distributions"
}

cat(justification, "\n\n")

# ========================================
# PART 2: BOOTSTRAP IMPLEMENTATION
# ========================================

cat("========================================\n")
cat("BOOTSTRAP UNCERTAINTY ESTIMATION\n")
cat("========================================\n\n")

set.seed(123)
B <- 10000  # Number of bootstrap samples

# Initialize storage for bootstrap estimates
boot_means <- numeric(B)
boot_vars <- numeric(B)
boot_medians <- numeric(B)
boot_q25 <- numeric(B)
boot_q75 <- numeric(B)
boot_cv <- numeric(B)  # Coefficient of variation

cat("Running", B, "bootstrap iterations...\n")

if(use_parametric){
  # Determine which parametric distribution to use
  if(shapiro_test$p.value > 0.05 && ks_normal$p.value > 0.05){
    cat("Using Normal distribution for parametric bootstrap\n\n")
    for(b in 1:B){
      boot_sample <- rnorm(length(type2_counts), mean = mean_val, sd = sd(type2_counts))
      boot_means[b] <- mean(boot_sample)
      boot_vars[b] <- var(boot_sample)
      boot_medians[b] <- median(boot_sample)
      boot_q25[b] <- quantile(boot_sample, 0.25)
      boot_q75[b] <- quantile(boot_sample, 0.75)
      boot_cv[b] <- sd(boot_sample) / mean(boot_sample)
    }
  } else {
    cat("Using Poisson distribution for parametric bootstrap\n\n")
    for(b in 1:B){
      boot_sample <- rpois(length(type2_counts), lambda = mean_val)
      boot_means[b] <- mean(boot_sample)
      boot_vars[b] <- var(boot_sample)
      boot_medians[b] <- median(boot_sample)
      boot_q25[b] <- quantile(boot_sample, 0.25)
      boot_q75[b] <- quantile(boot_sample, 0.75)
      boot_cv[b] <- sd(boot_sample) / mean(boot_sample)
    }
  }
} else {
  # Non-parametric bootstrap
  cat("Using Non-parametric (resampling) bootstrap\n\n")
  for(b in 1:B){
    boot_sample <- sample(type2_counts, size = length(type2_counts), replace = TRUE)
    boot_means[b] <- mean(boot_sample)
    boot_vars[b] <- var(boot_sample)
    boot_medians[b] <- median(boot_sample)
    boot_q25[b] <- quantile(boot_sample, 0.25)
    boot_q75[b] <- quantile(boot_sample, 0.75)
    boot_cv[b] <- sd(boot_sample) / mean(boot_sample)
  }
}

# Calculate bootstrap statistics
cat("BOOTSTRAP RESULTS (", B, " iterations)\n", sep="")
cat("----------------------------------------\n\n")

# Mean
cat("MEAN:\n")
cat("  Original estimate:", mean_val, "\n")
cat("  Bootstrap mean:", mean(boot_means), "\n")
cat("  Bootstrap SE:", sd(boot_means), "\n")
cat("  95% CI:", quantile(boot_means, c(0.025, 0.975)), "\n\n")

# Variance
cat("VARIANCE:\n")
cat("  Original estimate:", var_val, "\n")
cat("  Bootstrap mean:", mean(boot_vars), "\n")
cat("  Bootstrap SE:", sd(boot_vars), "\n")
cat("  95% CI:", quantile(boot_vars, c(0.025, 0.975)), "\n\n")

# Median
cat("MEDIAN:\n")
cat("  Original estimate:", median(type2_counts), "\n")
cat("  Bootstrap mean:", mean(boot_medians), "\n")
cat("  Bootstrap SE:", sd(boot_medians), "\n")
cat("  95% CI:", quantile(boot_medians, c(0.025, 0.975)), "\n\n")

# Coefficient of Variation
cat("COEFFICIENT OF VARIATION (SD/Mean):\n")
cat("  Original estimate:", sd(type2_counts)/mean_val, "\n")
cat("  Bootstrap mean:", mean(boot_cv), "\n")
cat("  Bootstrap SE:", sd(boot_cv), "\n")
cat("  95% CI:", quantile(boot_cv, c(0.025, 0.975)), "\n\n")

# ========================================
# PART 3: MONTE CARLO VALIDATION
# ========================================

cat("========================================\n")
cat("MONTE CARLO VALIDATION\n")
cat("========================================\n\n")

cat("Testing estimator robustness under different true distributions...\n\n")

set.seed(456)
n_mc <- 1000  # Number of Monte Carlo simulations
n_boot_mc <- 1000  # Bootstrap iterations per MC simulation
n_sample <- length(type2_counts)

# We'll test under 3 scenarios: Normal, Poisson, and Negative Binomial
scenarios <- list(
  Normal = list(name = "Normal", params = list(mean = mean_val, sd = sd(type2_counts))),
  Poisson = list(name = "Poisson", params = list(lambda = mean_val)),
  NegBinom = list(name = "Negative Binomial", params = list(size = 10, mu = mean_val))
)

mc_results <- list()

for(scenario_name in names(scenarios)){
  cat("Testing under", scenario_name, "distribution...\n")
  scenario <- scenarios[[scenario_name]]
  
  # Storage for MC results
  mc_mean_estimates <- numeric(n_mc)
  mc_mean_ses <- numeric(n_mc)
  mc_mean_coverage <- numeric(n_mc)
  
  mc_var_estimates <- numeric(n_mc)
  mc_var_ses <- numeric(n_mc)
  
  for(mc in 1:n_mc){
    # Generate data from true distribution
    if(scenario_name == "Normal"){
      true_data <- rnorm(n_sample, mean = scenario$params$mean, sd = scenario$params$sd)
    } else if(scenario_name == "Poisson"){
      true_data <- rpois(n_sample, lambda = scenario$params$lambda)
    } else if(scenario_name == "NegBinom"){
      true_data <- rnbinom(n_sample, size = scenario$params$size, mu = scenario$params$mu)
    }
    
    # Bootstrap on this MC sample (using non-parametric for robustness test)
    boot_means_mc <- numeric(n_boot_mc)
    boot_vars_mc <- numeric(n_boot_mc)
    
    for(b in 1:n_boot_mc){
      boot_sample <- sample(true_data, size = length(true_data), replace = TRUE)
      boot_means_mc[b] <- mean(boot_sample)
      boot_vars_mc[b] <- var(boot_sample)
    }
    
    # Store results
    mc_mean_estimates[mc] <- mean(true_data)
    mc_mean_ses[mc] <- sd(boot_means_mc)
    
    # Check coverage (does CI contain true parameter?)
    ci_mean <- quantile(boot_means_mc, c(0.025, 0.975))
    true_mean <- scenario$params[[ifelse(scenario_name == "Poisson", "lambda", "mean")]]
    mc_mean_coverage[mc] <- (ci_mean[1] <= true_mean) && (ci_mean[2] >= true_mean)
    
    mc_var_estimates[mc] <- var(true_data)
    mc_var_ses[mc] <- sd(boot_vars_mc)
  }
  
  mc_results[[scenario_name]] <- list(
    mean_estimates = mc_mean_estimates,
    mean_ses = mc_mean_ses,
    mean_coverage = mc_mean_coverage,
    var_estimates = mc_var_estimates,
    var_ses = mc_var_ses
  )
  
  cat("  Mean estimate bias:", mean(mc_mean_estimates - true_mean), "\n")
  cat("  Mean SE (avg):", mean(mc_mean_ses), "\n")
  cat("  95% CI coverage:", mean(mc_mean_coverage) * 100, "%\n")
  cat("  Variance estimate bias:", mean(mc_var_estimates - 
                                          ifelse(scenario_name == "Normal", scenario$params$sd^2,
                                                 ifelse(scenario_name == "Poisson", scenario$params$lambda,
                                                        scenario$params$mu + scenario$params$mu^2/scenario$params$size))), "\n\n")
}

# ========================================
# VISUALIZATION
# ========================================

cat("========================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("========================================\n\n")

# Create comprehensive visualization
pdf("bootstrap_analysis_results.pdf", width = 14, height = 10)

# Layout for bootstrap distributions
par(mfrow=c(2,3))

# 1. Bootstrap distribution of mean
hist(boot_means, breaks=50, prob=TRUE, main="Bootstrap Distribution of Mean",
     xlab="Mean", col="lightblue", border="white")
abline(v=mean_val, col="red", lwd=2, lty=2)
abline(v=quantile(boot_means, c(0.025, 0.975)), col="blue", lwd=2, lty=3)
legend("topright", legend=c("Original", "95% CI"), col=c("red", "blue"), lty=c(2,3), lwd=2)

# 2. Bootstrap distribution of variance
hist(boot_vars, breaks=50, prob=TRUE, main="Bootstrap Distribution of Variance",
     xlab="Variance", col="lightgreen", border="white")
abline(v=var_val, col="red", lwd=2, lty=2)
abline(v=quantile(boot_vars, c(0.025, 0.975)), col="blue", lwd=2, lty=3)

# 3. Bootstrap distribution of median
hist(boot_medians, breaks=50, prob=TRUE, main="Bootstrap Distribution of Median",
     xlab="Median", col="lightyellow", border="white")
abline(v=median(type2_counts), col="red", lwd=2, lty=2)
abline(v=quantile(boot_medians, c(0.025, 0.975)), col="blue", lwd=2, lty=3)

# 4. Bootstrap distribution of CV
hist(boot_cv, breaks=50, prob=TRUE, main="Bootstrap Distribution of CV",
     xlab="Coefficient of Variation", col="lightcoral", border="white")
abline(v=sd(type2_counts)/mean_val, col="red", lwd=2, lty=2)
abline(v=quantile(boot_cv, c(0.025, 0.975)), col="blue", lwd=2, lty=3)

# 5. Q-Q plot of bootstrap means
qqnorm(boot_means, main="Q-Q Plot: Bootstrap Means")
qqline(boot_means, col="red", lwd=2)

# 6. Bootstrap SE convergence
se_convergence <- sapply(seq(100, B, by=100), function(n) sd(boot_means[1:n]))
plot(seq(100, B, by=100), se_convergence, type="l", lwd=2,
     main="Bootstrap SE Convergence", xlab="Bootstrap Iterations", ylab="SE of Mean")
abline(h=sd(boot_means), col="red", lty=2)

# Monte Carlo validation plots
par(mfrow=c(2,3))

for(scenario_name in names(scenarios)){
  results <- mc_results[[scenario_name]]
  
  # Distribution of mean estimates
  hist(results$mean_estimates, breaks=30, prob=TRUE,
       main=paste("MC Mean Estimates\n", scenario_name),
       xlab="Mean Estimate", col="skyblue", border="white")
  
  # Distribution of SEs
  hist(results$mean_ses, breaks=30, prob=TRUE,
       main=paste("MC Standard Errors\n", scenario_name),
       xlab="Standard Error", col="lightpink", border="white")
}

# Coverage probability comparison
par(mfrow=c(1,1))
coverage_rates <- sapply(mc_results, function(x) mean(x$mean_coverage))
barplot(coverage_rates * 100, 
        main="95% Confidence Interval Coverage Rates\nAcross Different True Distributions",
        ylab="Coverage Rate (%)",
        ylim=c(0,100),
        col=c("steelblue", "coral", "seagreen"),
        names.arg=names(scenarios))
abline(h=95, col="red", lwd=2, lty=2)
text(x=1:3, y=coverage_rates*100 + 3, labels=paste0(round(coverage_rates*100, 1), "%"))

dev.off()

cat("Visualization saved to 'bootstrap_analysis_results.pdf'\n\n")

# ========================================
# SUMMARY AND RECOMMENDATIONS
# ========================================

cat("========================================\n")
cat("SUMMARY AND RECOMMENDATIONS\n")
cat("========================================\n\n")

cat("1. BOOTSTRAP METHOD:\n")
cat("  ", justification, "\n\n")

cat("2. UNCERTAINTY QUANTIFICATION:\n")
cat("   Mean: ", round(mean_val, 2), " ± ", round(sd(boot_means), 2), 
    " (95% CI: ", round(quantile(boot_means, 0.025), 2), "-", 
    round(quantile(boot_means, 0.975), 2), ")\n", sep="")
cat("   Variance: ", round(var_val, 2), " ± ", round(sd(boot_vars), 2),
    " (95% CI: ", round(quantile(boot_vars, 0.025), 2), "-", 
    round(quantile(boot_vars, 0.975), 2), ")\n", sep="")

cat("\n3. MONTE CARLO VALIDATION:\n")
cat("   The estimator shows robustness across different true distributions:\n")
for(scenario_name in names(scenarios)){
  cat("   - ", scenario_name, ": ", 
      round(mean(mc_results[[scenario_name]]$mean_coverage)*100, 1), 
      "% coverage\n", sep="")
}

optimal_coverage <- which.max(sapply(mc_results, function(x) mean(x$mean_coverage)))
cat("\n   Best performance under:", names(scenarios)[optimal_coverage], "\n")

if(min(coverage_rates) < 0.90){
  cat("\n   WARNING: Coverage below 90% for some distributions.\n")
  cat("   Consider increasing sample size or using different estimation method.\n")
} else {
  cat("\n   All coverage rates acceptable (>90%).\n")
  cat("   Bootstrap method provides reliable uncertainty estimates.\n")
}

cat("\n4. RECOMMENDATION:\n")
if(!use_parametric && all(coverage_rates > 0.93)){
  cat("   Non-parametric bootstrap is RECOMMENDED:\n")
  cat("   - Provides robust estimates across distributions\n")
  cat("   - No assumptions about data distribution needed\n")
  cat("   - Consistent performance in validation\n")
} else if(use_parametric && max(coverage_rates) > 0.95){
  cat("   Parametric bootstrap is RECOMMENDED:\n")
  cat("   - Data fits assumed distribution well\n")
  cat("   - More efficient (lower SE for same # iterations)\n")
  cat("   - Excellent coverage in validation\n")
} else {
  cat("   MIXED RESULTS - Consider:\n")
  cat("   - Collecting more data for better distribution identification\n")
  cat("   - Using non-parametric bootstrap as safer choice\n")
  cat("   - Consulting domain expert about data generation process\n")
}

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("========================================\n")
