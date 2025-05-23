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
