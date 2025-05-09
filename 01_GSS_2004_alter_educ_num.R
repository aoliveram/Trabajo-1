# Import dataset
gss_egor <- readRDS("trabajo_1_files/gss_egor.rds")

# Set seed for reproducibility
set.seed(1)

# Target frequencies from gss_alters$educ_num
target_freq <- c(14, 58, 194, 665, 611, 581, 393)
names(target_freq) <- c(4, 8, 11, 12, 14, 16, 18)
total_n <- sum(target_freq)

# Function to generate education years from a mixture of normal distributions
generate_educ_distribution <- function(n, weight1, mu1, sigma1, mu2, sigma2) {
  # Determine which component each sample comes from
  component <- runif(n) < weight1
  
  # Generate samples from each component
  samples <- numeric(n)
  samples[component] <- rnorm(sum(component), mu1, sigma1)
  samples[!component] <- rnorm(sum(!component), mu2, sigma2)
  
  # Round to nearest year and ensure within 0-20 range
  samples <- round(pmax(0, pmin(20, samples)))
  
  return(samples)
}

# Function to bin education years according to the mapping rules
bin_education_years <- function(years) {
  result <- numeric(length(years))
  
  # Apply the mapping rules based on the EDUC1 categories
  for (i in seq_along(years)) {
    year <- years[i]
    
    if (year >= 1 && year <= 6) {
      result[i] <- 4       # 1-6 years -> 4
    } else if (year >= 7 && year <= 9) {
      result[i] <- 8       # 7-9 years -> 8
    } else if (year >= 10 && year <= 11) {
      result[i] <- 11      # 10-11 years -> 11 (part of 10-12 no diploma)
    } else if (year == 12) {
      # For year 12, we need to distribute between "11" (no diploma) and "12" (HS grad)
      # Based on typical HS graduation rates, some people with 12 years 
      # of education might not have a diploma
      result[i] <- sample(c(11, 12), 1, prob = c(0.12, 0.88))
    } else if (year >= 13 && year <= 15) {
      result[i] <- 14      # 13-15 years -> 14 (Some college/Associate)
    } else if (year >= 16 && year <= 17) {
      result[i] <- 16      # 16-17 years -> 16 (Bachelor's)
    } else if (year >= 18 && year <= 20) {
      result[i] <- 18      # 18-20 years -> 18 (Graduate/professional)
    } else {
      # For values outside our expected range
      result[i] <- NA
    }
  }
  
  return(result)
}

# Parameters for the education distribution
# These parameters have been chosen to approximate the target frequencies
weight1 <- 0.65  # Weight for the first component (lower education)
mu1 <- 12.2      # Mean for the first component
sigma1 <- 2.8    # SD for the first component
mu2 <- 16.3      # Mean for the second component (higher education)
sigma2 <- 1.7    # SD for the second component

# Generate the synthetic education years
synthetic_years <- generate_educ_distribution(total_n, weight1, mu1, sigma1, mu2, sigma2)

hist(synthetic_years, breaks = 0:21, col = "lightblue",
     main = "Synthetic Education Distribution", 
     xlab = "Years of Education", 
     ylab = "Frequency")

# Bin the synthetic years according to the mapping rules
binned_values <- bin_education_years(synthetic_years)

hist(binned_values)

# Compare with target frequencies
synthetic_freq <- table(factor(binned_values, levels = names(target_freq)))
comparison <- rbind(target_freq, as.integer(synthetic_freq))
rownames(comparison) <- c("Target", "Synthetic")

print(comparison)

# Calculate chi-square goodness of fit
chi_sq <- sum((synthetic_freq - target_freq)^2 / target_freq)
print(paste("Chi-square goodness of fit:", round(chi_sq, 2)))

# --- Optimization -------------------------------------------------------------

# Objective function for optimization
objective <- function(params) {
  # Extract parameters
  weight1 <- params[1]
  mu1 <- params[2]
  sigma1 <- params[3]
  mu2 <- params[4]
  sigma2 <- params[5]
  
  # Generate synthetic distribution
  synthetic_years <- generate_educ_distribution(total_n, weight1, mu1, sigma1, mu2, sigma2)
  binned_values <- bin_education_years(synthetic_years)
  
  # Calculate error
  synthetic_freq <- table(factor(binned_values, levels = names(target_freq)))
  chi_sq <- sum((synthetic_freq - target_freq)^2 / target_freq)
  
  return(chi_sq)
}

# Initial parameters
init_params <- c(0.65, 12.2, 2.8, 16.3, 1.7)

# Parameter bounds
lower_bounds <- c(0.4, 10, 1, 14, 1)
upper_bounds <- c(0.9, 14, 4, 18, 3)

# Run optimization
opt_result <- optim(init_params, objective, 
                    method = "L-BFGS-B",
                    lower = lower_bounds,
                    upper = upper_bounds)

# Get optimized parameters
best_params <- opt_result$par
print(best_params)

# Generate distribution with best parameters
best_years <- generate_educ_distribution(total_n, best_params[1], best_params[2], 
                                         best_params[3], best_params[4], best_params[5])

# --- Comparing Distributions --------------------------------------------------

best_years_binned <- bin_education_years(best_years)

# Create overlapping histograms
# This would reproduce the style of plot shown in your question
par(mfrow=c(1,1))
hist(synthetic_years, breaks=0:21, col=rgb(1,0.4,0.4,0.5),
     main="Education Distribution Comparison",
     xlab="Years of Education",
     ylim=c(0, max(table(synthetic_years), table(gss_egos$educ_num))))
hist(best_years, breaks=0:21, col=rgb(0.2,0.4,0.6,0.5), add=TRUE)
abline(v=mean(synthetic_years), col=rgb(1,0.4,0.4,0.8), lty=2)
abline(v=mean(gss_egos$educ_num), col=rgb(0.2,0.4,0.6,0.8), lty=2)
legend("topright", c("Egos", "Alters (Synthetic)"), 
       fill=c(rgb(0.2,0.4,0.6,0.5), rgb(1,0.4,0.4,0.5)))


par(mfrow=c(1,1))
hist(best_years_binned, breaks=0:21, col=rgb(0.4,0.8,0.4,0.5))
hist(gss_alters$educ_num, breaks=0:21, col=rgb(0.2,0.4,0.6,0.5), add=TRUE)
# abline(v=mean(synthetic_years), col=rgb(1,0.4,0.4,0.8), lty=2)
# abline(v=mean(gss_egos$educ_num), col=rgb(0.2,0.4,0.6,0.8), lty=2)
 legend("topleft", c("best_years_binned", "gss_alters$educ_num"), 
        fill=c(rgb(0.4,0.8,0.4,0.5), rgb(0.2,0.4,0.6,0.5)))


