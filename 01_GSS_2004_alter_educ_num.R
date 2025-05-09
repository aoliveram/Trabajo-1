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


#-------------------------------------------------------------------------------
# PRUEBA CON 3 NORMALES 
#-------------------------------------------------------------------------------


# Función para generar educación desde una mezcla de TRES normales
generate_educ_distribution_3norm <- function(n, w1, w2, mu1, sd1, mu2, sd2, mu3, sd3) {
  # w3 se calcula para que la suma de los pesos sea 1
  w3 <- 1 - w1 - w2
  
  # Comprobaciones básicas de los pesos y sigmas
  if (w1 < 0 || w2 < 0 || w3 < 0 || w1 > 1 || w2 > 1 || w3 > 1 || (w1 + w2 > 1)) {
    # Si los pesos no son válidos, devuelve NAs o maneja el error
    # Esto es crucial para que `optim` no se rompa con parámetros inválidos
    # Devolver Inf en la función objetivo es mejor que NAs aquí.
    # Para la generación directa, podrías devolver NAs o parar con error.
    # Aquí, asumimos que los parámetros de entrada ya son válidos post-optimización.
  }
  if (sd1 <= 0 || sd2 <= 0 || sd3 <= 0) {
    # Sigmas deben ser positivas
  }
  
  component <- sample(1:3, size = n, replace = TRUE, prob = c(w1, w2, w3))
  
  samples <- numeric(n)
  samples[component == 1] <- rnorm(sum(component == 1), mu1, sd1)
  samples[component == 2] <- rnorm(sum(component == 2), mu2, sd2)
  samples[component == 3] <- rnorm(sum(component == 3), mu3, sd3)
  
  # Redondear y asegurar dentro del rango 0-20
  samples <- round(pmax(0, pmin(20, samples)))
  return(samples)
}

# Función objetivo para la optimización con 3 componentes normales
objective_3norm <- function(params) {
  w1 <- params[1]
  w2 <- params[2]
  mu1 <- params[3]
  sd1 <- params[4]
  mu2 <- params[5]
  sd2 <- params[6]
  mu3 <- params[7]
  sd3 <- params[8]
  
  # Penalizar parámetros inválidos para guiar al optimizador
  # Pesos deben estar entre 0 y 1, y w1+w2 <= 1
  if (w1 < 0 || w1 > 1 || w2 < 0 || w2 > 1 || (w1 + w2) > 1) return(Inf)
  # Sigmas deben ser positivas
  if (sd1 <= 0 || sd2 <= 0 || sd3 <= 0) return(Inf)
  # Medias dentro de un rango razonable (ej. 0-20)
  # if (mu1 < 0 || mu1 > 20 || mu2 < 0 || mu2 > 20 || mu3 < 0 || mu3 > 20) return(Inf)
  
  
  synthetic_years_3n <- generate_educ_distribution_3norm(total_n, w1, w2, mu1, sd1, mu2, sd2, mu3, sd3)
  
  # Si por alguna razón la generación falló (aunque las comprobaciones anteriores deberían evitarlo)
  if (any(is.na(synthetic_years_3n))) return(Inf) 
  
  # Usar tu función de binning consistente (¡asegúrate que es la correcta!)
  binned_values_3n <- bin_education_years(synthetic_years_3n) 
  
  synthetic_freq_3n <- table(factor(binned_values_3n, levels = names(target_freq)))
  
  # Asegurar que no haya NAs o longitudes incorrectas
  # (factor con levels ya debería manejar la mayoría de esto)
  current_target_freq <- target_freq
  # Si synthetic_freq_3n no tiene todos los niveles, table(factor(...)) lo arregla
  
  # Evitar división por cero en Chi-cuadrado
  # Considerar solo bins donde target_freq > 0
  valid_idx <- current_target_freq > 0
  
  chi_sq <- sum(((synthetic_freq_3n[valid_idx] - current_target_freq[valid_idx])^2) / current_target_freq[valid_idx])
  
  return(chi_sq)
}

# --- Configuración y Ejecución de la Optimización para 3 Normales ---
set.seed(124) # Para reproducibilidad en esta nueva optimización

# Parámetros iniciales para 3 componentes (w1, w2, mu1, sd1, mu2, sd2, mu3, sd3)
# Estos son solo ejemplos, ajústalos según tu intuición de los 3 picos:
# ej. HS, Some College/Assoc, Bach/Grad
init_params_3norm <- c(0.35, 0.30,  12.0, 1.5,  14.0, 1.0,  16.5, 1.8) 

# Límites para los parámetros
# (w1, w2 van de 0 a 1, pero su suma también debe ser <=1. L-BFGS-B maneja esto bien)
# (medias entre ~6 y ~19, sigmas > 0 y razonables)
# lower_bounds_3norm <- c(0.05, 0.05,  10.0, 0.5,  13.0, 0.5,  15.0, 0.5)
# upper_bounds_3norm <- c(0.70, 0.70,  13.0, 3.0,  15.0, 2.5,  19.0, 3.5)
lower_bounds_3norm <- c(0.05, 0.05,  10.5, 0.5,  13.0, 0.5,  15.5, 0.8)
upper_bounds_3norm <- c(0.60, 0.60,  12.8, 2.5,  15.2, 2.5,  19.5, 4.0)

cat("Iniciando optimización para 3 componentes normales...\n")
opt_result_3norm <- optim(
  par = init_params_3norm,
  fn = objective_3norm,
  method = "L-BFGS-B",
  lower = lower_bounds_3norm,
  upper = upper_bounds_3norm,
  control = list(maxit = 500, trace = 1) # maxit más alto, trace para ver progreso
)

cat("Optimización para 3 componentes completada.\n")
best_params_3norm <- opt_result_3norm$par
print("Mejores parámetros (3 normales):")
print(best_params_3norm)
print(paste("Valor del Chi-cuadrado (3 normales):", round(opt_result_3norm$value, 2)))

# Generar distribución final con los mejores parámetros de 3 normales
synthetic_years_best_3norm <- generate_educ_distribution_3norm(
  n = total_n,
  w1 = best_params_3norm[1], w2 = best_params_3norm[2],
  mu1 = best_params_3norm[3], sd1 = best_params_3norm[4],
  mu2 = best_params_3norm[5], sd2 = best_params_3norm[6],
  mu3 = best_params_3norm[7], sd3 = best_params_3norm[8]
)

synthetic_years_best_binned_3norm <- bin_education_years(synthetic_years_best_3norm)

# --- Comparación y Visualización para 3 Normales ---

# Comparación de frecuencias
comparison_3norm <- rbind(
  target_freq,
  as.integer(table(factor(synthetic_years_best_binned_3norm, levels = names(target_freq))))
)
rownames(comparison_3norm) <- c("Target (gss_alters$educ_num)", "Synthetic (3-norm binned)")
print("Comparación de frecuencias (3 normales):")
print(comparison_3norm)

# Histograma de comparación (similar al que ya tienes para el modelo de 2 normales)
# Primero, la distribución de los años sintéticos (0-20) contra los egos
par(mfrow=c(1,1)) # Asegurar una sola gráfica
hist(gss_egos$educ_num, breaks=seq(-0.5, 20.5, by=1), col=rgb(1,0.4,0.4,0.5), # Rojo para Egos
     main="Education Distribution Comparison (Egos vs. Synthetic 3-Norm)",
     xlab="Years of Education",
     ylab="Frequency",
     ylim=c(0, max(table(gss_egos$educ_num), table(synthetic_years_best_3norm))),
     xaxt="n") # Suprimir ejes x por defecto
axis(1, at=0:20, labels=0:20) # Eje x más limpio
hist(synthetic_years_best_3norm, breaks=seq(-0.5, 20.5, by=1), col=rgb(0.2,0.4,0.8,0.5), add=TRUE) # Azul para Alters sintéticos
abline(v=mean(gss_egos$educ_num, na.rm=TRUE), col=rgb(1,0.4,0.4,0.9), lty=2, lwd=2)
abline(v=mean(synthetic_years_best_3norm, na.rm=TRUE), col=rgb(0.2,0.4,0.8,0.9), lty=2, lwd=2)
legend("topleft", c("Egos (GSS)", "Alters (Synthetic 3-Norm)"),
       fill=c(rgb(1,0.4,0.4,0.5), rgb(0.2,0.4,0.8,0.5)))


# Segundo, la distribución binned sintética contra los alters originales binned
title_str_3norm <- sprintf(
  "Education Dist (Binned) Comparison (3-Norm)\n(w1=%.2f,w2=%.2f,w3=%.2f\nmu1=%.2f,s1=%.2f, mu2=%.2f,s2=%.2f, mu3=%.2f,s3=%.2f)",
  best_params_3norm[1], best_params_3norm[2], 1-best_params_3norm[1]-best_params_3norm[2],
  best_params_3norm[3], best_params_3norm[4], best_params_3norm[5], best_params_3norm[6],
  best_params_3norm[7], best_params_3norm[8]
)

# Calcular el ylim dinámicamente para el histograma binned
max_freq_binned_3norm <- max(c(target_freq, table(synthetic_years_best_binned_3norm)))

# Histograma de comparación para los datos binned
h1_data <- gss_alters$educ_num[!is.na(gss_alters$educ_num)]
h2_data <- synthetic_years_best_binned_3norm[!is.na(synthetic_years_best_binned_3norm)]

# Definir los breaks para que los bins sean los números exactos
unique_bins <- sort(unique(c(h1_data, h2_data)))
breaks_binned <- c(unique_bins - 0.5, max(unique_bins) + 0.5)

hist_gss_alters <- hist(h1_data, breaks = breaks_binned, plot = FALSE)
hist_synthetic_binned <- hist(h2_data, breaks = breaks_binned, plot = FALSE)

# Usar barplot para mejor control de superposición y etiquetas
barplot_data <- rbind(hist_synthetic_binned$counts, hist_gss_alters$counts)
# Puede que necesites ajustar el orden si quieres que gss_alters esté "debajo"
# barplot_data <- rbind(hist_gss_alters$counts, hist_synthetic_binned$counts)


# Crear el gráfico de barras superpuesto
bp <- barplot(barplot_data, beside = TRUE, # 'beside = FALSE' para apilado, 'TRUE' para lado a lado
              col = c(rgb(0.4,0.8,0.4,0.7), rgb(0.2,0.4,0.6,0.7)), # Verde para sintético, Azul para GSS
              main = title_str_3norm,
              xlab = "Binned Education Years (Alters)",
              ylab = "Frequency",
              ylim = c(0, max_freq_binned_3norm * 1.1),
              names.arg = hist_gss_alters$mids, # O usa los nombres de tus bins
              cex.names = 0.8)
legend("topright", c("Synthetic (3-Norm Binned)", "GSS Alters (Original Binned)"),
       fill=c(rgb(0.4,0.8,0.4,0.7), rgb(0.2,0.4,0.6,0.7)), cex=0.8)

# Si prefieres el estilo de histogramas superpuestos como antes:
# plot(hist_gss_alters, col=rgb(0.2,0.4,0.6,0.5), # Azul para GSS Alters
#      main=title_str_3norm,
#      xlab="Binned Education Years (Alters)",
#      ylab="Frequency",
#      ylim=c(0, max_freq_binned_3norm * 1.1),
#      xaxt = "n")
# axis(1, at = hist_gss_alters$mids, labels = names(target_freq)) # Usar los nombres de los bins como etiquetas
# plot(hist_synthetic_binned, col=rgb(0.4,0.8,0.4,0.5), add=TRUE) # Verde para sintético
# legend("topright", c("Synthetic (3-Norm Binned)", "GSS Alters (Original Binned)"),
#        fill=c(rgb(0.4,0.8,0.4,0.5), rgb(0.2,0.4,0.6,0.5)), cex=0.8)


cat("Script de 3 normales completado.\n")