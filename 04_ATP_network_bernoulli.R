library(dplyr)

ATP_W3_sub <- readRDS("trabajo_1_files/ATP_W3_imput.rds")

# sort(unique(ATP_W3_sub$age))
# unique(ATP_W3_sub$sex)
# sort(unique(ATP_W3_sub$educ_num))
# levels(ATP_W3_sub$race)
# levels(ATP_W3_sub$relig)

gss_egor <- readRDS("trabajo_1_files/gss_egor.rds")
gss_egos <- gss_egor$ego
gss_alters <- gss_egor$alter

# sort(unique(gss_egos$age))
# unique(gss_egos$sex)
# sort(unique(gss_egos$educ_num))
# levels(gss_egos$race)
# levels(gss_egos$relig)
# 
# sort(unique(gss_alters$age))
# unique(gss_alters$sex)
# sort(unique(gss_alters$educ_num))
# levels(gss_alters$race)
# levels(gss_alters$relig)

# ------------- P Matrix construction ------------------------------------------

N_atp <- nrow(ATP_W3_sub)

# --- Coeficientes del modelo de Smith et al. 2014 (Tabla 3, Modelo 2, All Ties) ---
# Estos son los Beta_0 a Beta_11 de la Ecuación 1 de McPherson & Smith (2019)
# Por ejemplo:
# beta_0_base = Intercepto del modelo GSS 1985
# beta_1_base = Coeficiente para Different_Race en 1985
# ...
# beta_6 = Coeficiente para el término Year (Year=1 para 2004)
# beta_7 = Coeficiente para la interacción Different_Race * Year
# ...

beta_coefficients <- c(
  Intercept          = -14.519,
  Different_Race     = -1.959,
  Different_Religion = -1.270,
  Different_Sex      = -0.373,
  Age_Difference     = -0.047,
  Education_Difference = -0.157,
  Year_2004          = -0.052, # beta_6
  Diff_Race_x_Year   =  0.264, # beta_7
  Diff_Relig_x_Year  = -0.215, # beta_8
  Diff_Sex_x_Year    =  0.144, # beta_9
  Age_Diff_x_Year    = -0.005, # beta_10
  Educ_Diff_x_Year   = -0.044  # beta_11
)

# Convertimos ATP_W3_sub$sex, que es un factor "Male"/"Female", a numeric:
if (is.factor(ATP_W3_sub$sex)) {
  print('sex is factor. Converting to numeric ATP_W3_sub$sex_numeric : Male=1, Female=2')
  ATP_W3_sub$sex_numeric <- as.numeric(ATP_W3_sub$sex) # Male=1, Female=2 si Male es el primer nivel
} else {
  ATP_W3_sub$sex_numeric <- ATP_W3_sub$sex
}

# Una matriz P completa de 4000x4000 es 4000*4000*8 bytes (double) ~ 128MB, manejable.
P_matrix <- matrix(0.0, nrow = N_atp, ncol = N_atp)
rownames(P_matrix) <- ATP_W3_sub$QKEY # Asumiendo que tienes un ID único
colnames(P_matrix) <- ATP_W3_sub$QKEY

# --- Calcular p_ij para cada par (i, j) ---

# Asignar Year_ij = 1 (refleja el estado más reciente, valores GSS 2004)
Year_ij_val <- 1

# Coeficientes efectivos para Year_ij = 1
beta_eff_intercept <- beta_coefficients["Intercept"] + beta_coefficients["Year_2004"] * Year_ij_val
beta_eff_race    <- beta_coefficients["Different_Race"] + beta_coefficients["Diff_Race_x_Year"] * Year_ij_val
beta_eff_relig   <- beta_coefficients["Different_Religion"] + beta_coefficients["Diff_Relig_x_Year"] * Year_ij_val
beta_eff_sex     <- beta_coefficients["Different_Sex"] + beta_coefficients["Diff_Sex_x_Year"] * Year_ij_val
beta_eff_age     <- beta_coefficients["Age_Difference"] + beta_coefficients["Age_Diff_x_Year"] * Year_ij_val
beta_eff_educ    <- beta_coefficients["Education_Difference"] + beta_coefficients["Educ_Diff_x_Year"] * Year_ij_val

# Calculando la matriz P de probabilidades de lazo...
for (i in 1:N_atp) {
  if (i %% 100 == 0) cat("Procesando fila i =", i, "de", N_atp, "\n")
  for (j in (i + 1):N_atp) { # Solo calcular para la mitad superior (j > i)
    if (j > N_atp) next 
    
    # Calcular diferencias demográficas para el par (i, j)
    RaceDiff_ij <- ifelse(ATP_W3_sub$race[i] != ATP_W3_sub$race[j], 1, 0)
    ReligionDiff_ij <- ifelse(ATP_W3_sub$relig[i] != ATP_W3_sub$relig[j], 1, 0)
    SexDiff_ij <- ifelse(ATP_W3_sub$sex_numeric[i] != ATP_W3_sub$sex_numeric[j], 1, 0)
    AgeDiff_ij <- abs(ATP_W3_sub$age[i] - ATP_W3_sub$age[j])
    EducDiff_ij <- abs(ATP_W3_sub$educ_num[i] - ATP_W3_sub$educ_num[j])
    
    # Calcular eta_ij y p_ij
    if (any(is.na(c(RaceDiff_ij, ReligionDiff_ij, SexDiff_ij, AgeDiff_ij, EducDiff_ij)))) {
      # un 4% de ellos son NA.
      # > sum(is.na(P_matrix))/sum(!is.na(P_matrix))
      #[1] 0.04790112
      p_ij <- 0
    } else {
      # Calcular eta_ij
      eta_ij <- beta_eff_intercept +
        beta_eff_race    * RaceDiff_ij +
        beta_eff_relig   * ReligionDiff_ij +
        beta_eff_sex     * SexDiff_ij +
        beta_eff_age     * AgeDiff_ij +
        beta_eff_educ    * EducDiff_ij
      # calcular p_ij
      p_ij <- 1 / (1 + exp(-eta_ij))
    }
    
    P_matrix[i, j] <- p_ij
    P_matrix[j, i] <- p_ij # Matriz simétrica
  }
}

diag(P_matrix) <- 0 # No hay lazos de un nodo consigo mismo

# Guardamos matriz
# saveRDS(P_matrix, 'trabajo_1_files/P_Matrix.rds')
