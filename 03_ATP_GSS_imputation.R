library(dplyr)

# Cargamos datos GSS 2004
gss_egor <- readRDS("trabajo_1_files/gss_egor.rds")
gss_egos <- gss_egor$ego
gss_alters <- gss_egor$alter

hist(gss_egos$age)
hist(gss_egos$educ_num)

# Cargamos datos ATP W3

ATP_W3_sub <- readRDS("trabajo_1_files/ATP_W3_sub.rds")

hist(ATP_W3_sub$F_AGECAT_TYPOLOGY) # Edad

# Qué datos tenemos ---------------------------------

# Edad
labels[["F_AGECAT_TYPOLOGY"]]
ATP_W3_sub$F_AGECAT_TYPOLOGY

# Educación
labels[["F_EDUCCAT_TYPOLOGY"]]
ATP_W3_sub$F_EDUCCAT_TYPOLOGY

# Race
labels[["F_RACETHN_TYPOLOGY"]]
ATP_W3_sub$F_RACETHN_TYPOLOGY

# Sex
labels[["F_SEX_FINAL"]]
ATP_W3_sub$F_SEX_FINAL

# Religion
labels[["F_RELIG_TYPOLOGY"]]
ATP_W3_sub$F_RELIG_TYPOLOGY

# Ver las primeras filas del data frame
head(ATP_W3)
head(ATP_W3_sub)

# ------------------------------------------------------------------------------
# EDAD
# ------------------------------------------------------------------------------

# --- Veamos cómo se comparan las EDADES en GSS y en ATP -----------------------

# ]17,29], ]29,49], ]49,64], ]65, inf]
age_breaks_gss <- c(17, 29, 49, 64, Inf) 
age_labels_gss <- c("18-29", "30-49", "50-64", "65+")

# Proporciones de estas categorías en GSS
gss_egos$age_category_gss <- cut(gss_egos$age,
                                 breaks = age_breaks_gss,
                                 labels = age_labels_gss,
                                 right = TRUE, # Intervalos son (lim_inf, lim_sup]
                                 include.lowest = TRUE)

gss_age_proportions <- prop.table(table(gss_egos$age_category_gss))
print(gss_age_proportions)

# Proporciones de las categorías de edad en ATP  -(Ya es un factor con etiquetas de haven)
ATP_W3_sub$F_AGECAT_TYPOLOGY_factor <- haven::as_factor(ATP_W3$F_AGECAT_TYPOLOGY[ATP_W3$QKEY %in% ATP_W3_sub$QKEY])

atp_age_proportions <- prop.table(table(ATP_W3_sub$F_AGECAT_TYPOLOGY_factor))
print(atp_age_proportions)

# Comparación VISUAL

ages_gss_atp <- rbind(
  data.frame(source = "GSS", category = names(gss_age_proportions), proportion = as.numeric(gss_age_proportions)),
  data.frame(source = "ATP", category = names(atp_age_proportions), proportion = as.numeric(atp_age_proportions))
)

print(ggplot(ages_gss_atp, aes(x = category, y = proportion, fill = source)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Comparación de Distribuciones de Categorías de Edad (GSS vs ATP)",
             x = "Categoría de Edad", y = "Proporción") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Comparación Chi-Cuadrado
# Necesitamos las cuentas, no las proporciones
gss_counts <- table(gss_egos$age_category_gss)
atp_counts_raw <- table(ATP_W3_sub$F_AGECAT_TYPOLOGY_factor)

# tabla de contingencia para el test
contingency_table_age <- rbind(GSS = as.numeric(gss_counts_aligned), ATP = as.numeric(atp_counts_aligned))
colnames(contingency_table_age) <- all_age_cats

print(contingency_table_age)

# Test Chi-cuadrado
chi_sq_test_age <- chisq.test(contingency_table_age)
print(chi_sq_test_age)

if (chi_sq_test_age$p.value < 0.05) {
  cat("Conclusión: Las distribuciones de categorías de edad entre GSS y ATP son significativamente DIFERENTES (p < 0.05).\n")
  cat("Sin embargo, el muestreo *dentro* de los estratos definidos por las categorías ATP sigue siendo una estrategia razonable.\n")
} else {
  cat("Conclusión: No hay evidencia de una diferencia significativa en las distribuciones de categorías de edad (p >= 0.05).\n")
  cat("Esto apoya la idea de que la distribución de edad categórica de GSS es comparable a la de ATP.\n")
}

# ---------- Imputación de edades en ATP ---------------------------------------

# PREPARACIÓN para la imputación
# Esto crea una lista donde cada elemento es un vector de edades de GSS para esa categoría.
gss_ages_by_category <- tapply(gss_egos$age, gss_egos$age_category_gss, list)

length(gss_ages_by_category$`18-29`)+length(gss_ages_by_category$`30-49`)+length(gss_ages_by_category$`50-64`)+length(gss_ages_by_category$`65+`)
nrow(gss_egos)

print(lapply(gss_ages_by_category, head, 5))

# Crear la nueva columna en ATP_W3_sub
ATP_W3_sub$age <- NA_integer_ # Inicializar con NA

# Establecer una semilla para reproducibilidad de la imputación
set.seed(123)

# Iterar para cada individuo de ATP
for (i in 1:nrow(ATP_W3_sub)) {
  # Obtenmos su categoría de edad
  atp_category <- as.character(ATP_W3_sub$F_AGECAT_TYPOLOGY_factor[i])
  
  # Verificar si la categoría existe
  if (atp_category %in% names(gss_ages_by_category)) {
    # vector de edades GSS para esa categoría ---> Aquí ya está la distribución de GSS !
    possible_gss_ages <- gss_ages_by_category[[atp_category]]
    
    imputed_age <- sample(possible_gss_ages, 1)
    ATP_W3_sub$age[i] <- imputed_age
    
  } else {
    warning(paste("Categoría de edad ATP no mapeada:", atp_category, "en la fila", i))
    ATP_W3_sub$age[i] <- NA # O alguna otra estrategia
  }
}

summary(ATP_W3_sub$age)
hist(ATP_W3_sub$age, main = "Edades Imputadas en ATP \n(basado en distribución GSS)", xlab = "Edad Imputada")


# ------------------------------------------------------------------------------
# EDUCACIÓN
# ------------------------------------------------------------------------------

# --- Veamos cómo se comparan los niveles de EDUCACIÖN en GSS y en ATP ---------

# Crear categoría colapsada en GSS
gss_egos$degree_atp3cat <- case_when(
  gss_egos$degree %in% c("Bach", "Grad") ~ "College graduate+",
  gss_egos$degree == "Assoc" ~ "Some college",
  gss_egos$degree %in% c("HS", "LTHS") ~ "HS graduate or less",
  TRUE ~ NA_character_
)
gss_egos$degree_atp3cat <- factor(gss_egos$degree_atp3cat,
                                  levels = c("College graduate+", "Some college", "HS graduate or less"))

# ATP ya tiene la variable F_EDUCCAT_TYPOLOGY: 1=College graduate+, 2=Some college, 3=HS graduate or less
ATP_W3_sub$F_EDUCCAT_TYPOLOGY_factor <- factor(
  ATP_W3_sub$F_EDUCCAT_TYPOLOGY,
  levels = c(1, 2, 3),
  labels = c("College graduate+", "Some college", "HS graduate or less")
)

# Tabla de contingencia
gss_counts <- table(gss_egos$degree_atp3cat)
atp_counts <- table(ATP_W3_sub$F_EDUCCAT_TYPOLOGY_factor)
contingency_table_educ <- rbind(GSS = as.numeric(gss_counts), ATP = as.numeric(atp_counts))
colnames(contingency_table_educ) <- levels(gss_egos$degree_atp3cat)
print(contingency_table_educ)

# Test Chi-cuadrado
chi_sq_test_educ <- chisq.test(contingency_table_educ)
print(chi_sq_test_educ)

if (chi_sq_test_educ$p.value < 0.05) {
  cat("Conclusión: Las distribuciones de categorías de educación entre GSS y ATP son significativamente DIFERENTES (p < 0.05).\n")
  cat("Sin embargo, el muestreo *dentro* de los estratos definidos por las categorías ATP sigue siendo una estrategia razonable.\n")
} else {
  cat("Conclusión: No hay evidencia de una diferencia significativa en las distribuciones de categorías de educación (p >= 0.05).\n")
  cat("Esto apoya la idea de que la distribución de educación categórica de GSS es comparable a la de ATP.\n")
}

# ---------- Imputación de educ_num en ATP -------------------------------------

# Creando una lista por cada categoría colapsada
gss_educ_by_cat <- split(gss_egos$educ_num, gss_egos$degree_atp3cat)

# Inicializamos columna 
ATP_W3_sub$educ_num <- NA_real_

set.seed(123)
for (i in 1:nrow(ATP_W3_sub)) {
  atp_cat <- as.character(ATP_W3_sub$F_EDUCCAT_TYPOLOGY_factor[i])
  if (atp_cat %in% names(gss_educ_by_cat)) {
    possible_vals <- gss_educ_by_cat[[atp_cat]]
    ATP_W3_sub$educ_num[i] <- sample(possible_vals, 1)
  } else {
    ATP_W3_sub$educ_num[i] <- NA
  }
}

summary(ATP_W3_sub$educ_num)
hist(ATP_W3_sub$educ_num, main = "Años de educación imputados en ATP", xlab = "Años de educación")
