library(ergm) # Para simulate.ergm y network
library(dplyr)

# --- 0. Preparación ---

ATP_W3_sub <- readRDS("trabajo_1_files/ATP_W3_imput.rds")
N_atp <- nrow(ATP_W3_sub)

gss_egor <- readRDS("trabajo_1_files/gss_egor.rds")
gss_egos <- gss_egor$ego
gss_alters <- gss_egor$alter

# Nos aseguramos que las categorías sean factores
if (!is.factor(ATP_W3_sub$sex)) ATP_W3_sub$sex <- factor(ATP_W3_sub$sex, labels = c("Male", "Female")) # Ajusta si es necesario
if (!is.factor(ATP_W3_sub$race)) ATP_W3_sub$race <- factor(ATP_W3_sub$race)
if (!is.factor(ATP_W3_sub$relig)) ATP_W3_sub$relig <- factor(ATP_W3_sub$relig)

# Creamos el objeto NEtwork (sin lazos, solo nodos y atributos)
atp_base_network <- network.initialize(N_atp, directed = FALSE)
set.vertex.attribute(atp_base_network, "age", ATP_W3_sub$age)
set.vertex.attribute(atp_base_network, "sex", as.character(ATP_W3_sub$sex)) # ergm a veces prefiere character para nodematch
set.vertex.attribute(atp_base_network, "educ_num", ATP_W3_sub$educ_num)
set.vertex.attribute(atp_base_network, "race", as.character(ATP_W3_sub$race))
set.vertex.attribute(atp_base_network, "relig", as.character(ATP_W3_sub$relig))

# Fórmula ERGM (el término 'edge' se iterará)
formula_homofilia_only <- ~ nodematch("race") +
                            nodematch("sex") +
                            absdiff("age") +
                            absdiff("educ_num") +
                            nodematch("relig")
# Si se añade gwesp, también necesitaría un coeficiente fijo

# Smith et al. modelan P(lazo | Diferencia), no P(lazo | Coincidencia).
# El signo de los coeficientes en ergm para 'nodematch' será POSITIVO si la similitud aumenta la probabilidad.
# Los coeficientes de Smith et al. son para 'Different_X', entonces para 'nodematch.X' el signo se invierte.

# Coeficientes de Smith (2014), Tabla 3, Modelo 2 "All Ties" - AJUSTADOS PARA 'nodematch'
# Year_ij_val = 1 (para ATP, análogo a 2004 GSS)
beta_s2014_raw <- c(Different_Race = -1.959, Different_Religion = -1.270, Different_Sex = -0.373,
                    Age_Difference = -0.047, Education_Difference = -0.157, Year_2004 = -0.052,
                    Diff_Race_x_Year = 0.264, Diff_Relig_x_Year = -0.215, Diff_Sex_x_Year = 0.144,
                    Age_Diff_x_Year = -0.005, Educ_Diff_x_Year = -0.044) # AGREGAR TAMBIÉN AQUÍ 'YEAR' ??

# Coeficientes efectivos para ATP (asumimos que 2004 es lo más cercano)
coef_eff_diff_race  <- beta_s2014_raw["Different_Race"] + beta_s2014_raw["Diff_Race_x_Year"]
coef_eff_diff_relig <- beta_s2014_raw["Different_Religion"] + beta_s2014_raw["Diff_Relig_x_Year"]
coef_eff_diff_sex   <- beta_s2014_raw["Different_Sex"] + beta_s2014_raw["Diff_Sex_x_Year"]
coef_eff_absdiff_age  <- beta_s2014_raw["Age_Difference"] + beta_s2014_raw["Age_Diff_x_Year"]
coef_eff_absdiff_educ <- beta_s2014_raw["Education_Difference"] + beta_s2014_raw["Educ_Diff_x_Year"]

# El coef de Smith et al. (-1.695 para raza efectiva en 2004) es la log-odds de lazar
# si eres de DIFERENTE raza comparado con la misma raza (si el intercepto captura la misma raza).
# Entonces, para nodematch("race") en ergm, el coeficiente debería ser POSITIVO si la misma raza es más probable.
# Es decir, -1 * (coef_eff_diff_race)
# Los términos `absdiff` ya tienen el signo correcto (negativo indica que mayor diferencia disminuye lazos)

coef_homofilia_fijos_ergm <- c(
  -coef_eff_diff_race,  # Para nodematch("race")
  -coef_eff_diff_sex,   # Para nodematch("sex")
  coef_eff_absdiff_age, # Para absdiff("age")
  coef_eff_absdiff_educ,# Para absdiff("educ_num")
  -coef_eff_diff_relig  # Para nodematch("relig")
)

# Asignamos nombres para ergm (REVISAR EN SIMULATE)
names(coef_homofilia_fijos_ergm) <- c("nodematch.race", "nodematch.sex",
                                      "absdiff.age", "absdiff.educ_num", "nodematch.relig")

# --- 1. Establecer Parámetros de Búsqueda ---
target_density <- 0.029
lower_bound_edges <- -10.0  # Límite inferior para el coeficiente de edges
upper_bound_edges <- -2.0   # Límite superior

max_iterations <- 20       # Número máximo de iteraciones
density_tolerance <- 0.001 # tolerancia al target_density
num_sim_per_iteration <- 5 # Número de redes a simular en cada paso para promediar densidad

control_sim_formula <- control.simulate.formula( # REVISAR convergencia
  MCMC.burnin = 1000 * N_atp,
  MCMC.interval = 100 * N_atp
)

# --- 2. Bucle de pptimización para 'edges' ---

# Variables de atributos para ergm
attribute_vars <- c("race", "sex", "age", "educ_num", "relig")
sapply(ATP_W3_sub[, attribute_vars], function(col) sum(is.na(col))) # verificar NA's
sapply(ATP_W3_sub[, attribute_vars], function(col) mean(is.na(col)) * 100) # porcentaje NA's
#ATP_W3_sub_2 <- ATP_W3_sub
ATP_W3_sub <- ATP_W3_sub[complete.cases(ATP_W3_sub[, attribute_vars]), ] # Quitamos NA's
N_atp <- nrow(ATP_W3_sub) # Número de individuos en ATP

# Cremos base Network (sin lazos, solo información de nodos y atributos)
atp_base_network <- network.initialize(N_atp, directed = FALSE)
set.vertex.attribute(atp_base_network, "age", ATP_W3_sub$age)
set.vertex.attribute(atp_base_network, "sex", as.character(ATP_W3_sub$sex))
set.vertex.attribute(atp_base_network, "educ_num", ATP_W3_sub$educ_num)
set.vertex.attribute(atp_base_network, "race", as.character(ATP_W3_sub$race))
set.vertex.attribute(atp_base_network, "relig", as.character(ATP_W3_sub$relig))

# Agregamos variables estructurales a la fórmula
full_formula_ergm <- update.formula(formula_homofilia_only, paste("~ edges + ."))
# Objeto para alojar coeficientes calibrados
calibrated_coef_edges <- NA

# Iteración para ajustar 'edges' a densidad objetivo target_density
for (iter in 1:max_iterations) {
  current_coef_edges <- (lower_bound_edges + upper_bound_edges) / 2
  
  current_coefs_named <- c(edges = current_coef_edges, coef_homofilia_fijos_ergm)
  
  cat(paste("Iteración", iter, ": Probando coef_edges =", round(current_coef_edges, 4), "\n"))
  
  sim_networks_list <- tryCatch({
    simulate(
      full_formula_ergm,
      basis = atp_base_network, # El objeto network base
      nsim = num_sim_per_iteration,
      coef = current_coefs_named,
      control = control_sim_formula,
      verbose = FALSE
    )
  }, error = function(e) {
    cat("Error durante la simulación en la iteración", iter, ":\n", conditionMessage(e), "\n")
    cat("Coeficientes actuales:", paste(names(current_coefs_named), round(current_coefs_named,3), collapse=", "), "\n")
    cat("Nombres esperados por la fórmula (aproximado):",paste(attr(terms(full_formula_ergm),"term.labels"),collapse=", "),"\n")
    return(NULL)
  })
  
  # Si falla la simulación
  if (is.null(sim_networks_list)) {
    cat("Fallo en la simulación, ajustando límites de búsqueda o deteniendo.\n")
    # Lógica simple para evitar quedarse atascado si el modelo es inestable con ciertos coefs
    if (abs(current_coef_edges - lower_bound_edges) < 1e-6 || abs(current_coef_edges - upper_bound_edges) < 1e-6) {
      warning("Los límites de búsqueda están muy juntos y la simulación sigue fallando.")
      break 
    }
    # Si current_coef_edges es más positivo, tendía a dar más lazos. Densidad es "mala" (ej. demasiado alta o el modelo es degenerado)
    # Si falló, quizás es demasiado positivo. Intentamos hacerlo más negativo.
    upper_bound_edges <- current_coef_edges 
    next
  }
  
  # Calculamos densidad media de las redes
  sim_densities <- sapply(sim_networks_list, network.density)
  avg_sim_density <- mean(sim_densities, na.rm = TRUE)
  
  cat(paste("   Densidad promedio simulada:", round(avg_sim_density, 5), "\n"))
  
  # Comparamos con density_tolerance
  if (abs(avg_sim_density - target_density) < density_tolerance) {
    calibrated_coef_edges <- current_coef_edges
    cat(paste("Convergencia alcanzada en iteración", iter, "! Coef_edges calibrado =", round(calibrated_coef_edges, 4), "\n"))
    break
  }
  
  # Ajustamos nuevos bordes para parámetros 
  if (avg_sim_density < target_density) {
    lower_bound_edges <- current_coef_edges
  } else {
    upper_bound_edges <- current_coef_edges
  }
  
  # máximo iteraciones
  if (iter == max_iterations) {
    calibrated_coef_edges <- current_coef_edges
    cat(paste("Máximo de iteraciones alcanzado. Coef_edges final (aproximado) =", round(calibrated_coef_edges, 4), "\n"))
    cat(paste("   Densidad obtenida con este coeficiente:", round(avg_sim_density, 5), "\n"))
  }
}

if (!is.na(calibrated_coef_edges)) {
  cat("\nCalibración exitosa para 'edges'.\n")
  final_ergm_coefs <- c(edges = calibrated_coef_edges, coef_homofilia_fijos_ergm)
  print(final_ergm_coefs)
} else {
  cat("\nNo se pudo calibrar el coeficiente 'edges'.\n")
}

# --- 3. Simulación ---

library(network)
library(sna)

# Ejemplo de simulación de UNA red final:
ATP_network_final_simulated <- simulate(
  full_formula_ergm,
  basis = atp_base_network, # --> base con 3168 nodos
  nsim = 1,
  coef = final_ergm_coefs,
  control = control_sim_formula,
  verbose = TRUE
)

# Estadísticos Básicos
print(summary(ATP_network_final_simulated))

network.size(ATP_network_final_simulated)
network.edgecount(ATP_network_final_simulated)
network.density(ATP_network_final_simulated)

grados <- degree(ATP_network_final_simulated, gmode="graph") # 'graph' para no dirigida
print(summary(grados))
hist(grados, main="Distribución de Grado de la Red Simulada", xlab="Grado", ylab="Frecuencia", breaks=50)
