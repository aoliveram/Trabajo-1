# Import dataset
gss_egor <- readRDS("trabajo_1_files/gss_egor.rds")
gss_egos <- gss_egor$ego
gss_alters <- gss_egor$alter

table(gss_egos$degree)
table(gss_egos$educ_num)

table(gss_alters$educ_cat)

# Set seed for reproducibility
set.seed(1)

# Target frequencies from table(gss_alters$educ_num)
target_freq <- c(14, 58, 194, 665, 611, 581, 393)
names(target_freq) <- c(4, 8, 11, 12, 14, 16, 18)
total_n <- sum(target_freq)

# Function to generate education years from a mixture of two normal distributions
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
# bin_education_years <- function(years) {
#   result <- numeric(length(years))
#   
#   # Apply the mapping rules based on the EDUC1 categories
#   for (i in seq_along(years)) {
#     year <- years[i]
#     
#     if (year >= 1 && year <= 6) {
#       result[i] <- 4       # 1-6 years -> 4
#     } else if (year >= 7 && year <= 9) {
#       result[i] <- 8       # 7-9 years -> 8
#     } else if (year >= 10 && year <= 11) {
#       result[i] <- 11      # 10-11 years -> 11 (part of 10-12 no diploma)
#     } else if (year == 12) {
#       # For year 12, we need to distribute between "11" (no diploma) and "12" (HS grad)
#       # Based on typical HS graduation rates, some people with 12 years 
#       # of education might not have a diploma
#       result[i] <- sample(c(11, 12), 1, prob = c(0.12, 0.88))
#     } else if (year >= 13 && year <= 15) {
#       result[i] <- 14      # 13-15 years -> 14 (Some college/Associate)
#     } else if (year >= 16 && year <= 17) {
#       result[i] <- 16      # 16-17 years -> 16 (Bachelor's)
#     } else if (year >= 18 && year <= 20) {
#       result[i] <- 18      # 18-20 years -> 18 (Graduate/professional)
#     } else {
#       # For values outside our expected range
#       result[i] <- NA
#     }
#   }
#   
#   return(result)
# }

# bin_education_years_simple_match
bin_education_years <- function(years) {
  sapply(years, function(year) {
    if (is.na(year)) return(NA)
    if (year >= 0 && year <= 6) return(4)    # Midpoint approx for 1-6 (cat 0)
    if (year >= 7 && year <= 9) return(8)    # Midpoint approx for 7-9 (cat 1)
    if (year >= 10 && year <= 11) return(11) # Represents 10-12 years (cat 2) *sin ser HS grad*
                                             # Si educ_cat 2 (10-12 años) realmente incluía a algunos con 12 años sin diploma,
                                             # esta regla es la que se usó para gss_alters$educ_num -> 11.
    if (year == 12) return(12)               # HS Grad (cat 3)
    if (year >= 13 && year <= 15) return(14) # Some College (cat 4) or Assoc (cat 5)
    if (year >= 16 && year <= 17) return(16) # Bach (cat 6)
    if (year >= 18 && year <= 20) return(18) # Grad (cat 7)
    return(NA) # Catch-all for unexpected values
  })
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
  
  # ... (generación de synthetic_years y binned_values) ...
  synthetic_freq <- table(factor(binned_values, levels = names(target_freq)))
  
  # Asegurarse de que target_freq y synthetic_freq tengan el mismo orden y longitud
  # y no haya ceros en target_freq para la división
  valid_bins <- names(target_freq)[target_freq > 0]
  current_target_freq <- target_freq[valid_bins]
  current_synthetic_freq <- synthetic_freq[valid_bins]
  
  # Pesos (ejemplo: dar más peso al bin "12")
  weights <- rep(1, length(current_target_freq))
  if ("12" %in% names(current_target_freq)) {
    weights[names(current_target_freq) == "12"] <- 2 # Pondera el error del bin 12 el doble
  }
  
  chi_sq_contributions <- ((current_synthetic_freq - current_target_freq)^2 / current_target_freq)
  weighted_chi_sq <- sum(weights * chi_sq_contributions)
  
  return(weighted_chi_sq)
  
  # # Calculate error
  # synthetic_freq <- table(factor(binned_values, levels = names(target_freq)))
  # chi_sq <- sum((synthetic_freq - target_freq)^2 / target_freq)
  # 
  # return(chi_sq)
}

# Initial parameters (weight1, mu1, sigma1, mu2, sigma2)
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
synthetic_years_best <- generate_educ_distribution(total_n, best_params[1], best_params[2], 
                                         best_params[3], best_params[4], best_params[5])

# --- Comparing Distributions --------------------------------------------------

synthetic_years_best_binned <- bin_education_years(synthetic_years_best)

# Create overlapping histograms
# This would reproduce the style of plot shown in your question
par(mfrow=c(1,1))
hist(gss_egos$educ_num, breaks=0:21, col=rgb(1,0.4,0.4,0.5),
     main="Education Distribution Comparison",
     xlab="Years of Education",
     ylim=c(0, max(table(synthetic_years_best), table(gss_egos$educ_num))))
hist(synthetic_years_best, breaks=0:21, col=rgb(0.2,0.4,0.6,0.5), add=TRUE)
abline(v=mean(gss_egos$educ_num), col=rgb(0.2,0.4,0.6,0.8), lty=2)
abline(v=mean(synthetic_years_best), col=rgb(1,0.4,0.4,0.8), lty=2)
legend("topleft", c("Egos", "Alters (Synthetic)"), 
       fill=c(rgb(1,0.4,0.4,0.5), rgb(0.2,0.4,0.6,0.5)))


par(mfrow=c(1,1))
hist(synthetic_years_best_binned, breaks=0:21, col=rgb(0.4,0.8,0.4,0.5),
     main = paste("Education Distribution (binned) Comparison\n(w1:", 
                  sprintf("%.2f", best_params[1]), ", mu1:", sprintf("%.2f", best_params[2]),
                  ", sigma1:", sprintf("%.2f", best_params[3]), ", mu2:", sprintf("%.2f", best_params[4]),
                  ", sigma2:", sprintf("%.2f", best_params[5]), ")"))
hist(gss_alters$educ_num, breaks=0:21, col=rgb(0.2,0.4,0.6,0.5), add=TRUE)
legend("topleft", c("synthetic_years_best_binned", "gss_alters$educ_num"), 
        fill=c(rgb(0.4,0.8,0.4,0.5), rgb(0.2,0.4,0.6,0.5)))


