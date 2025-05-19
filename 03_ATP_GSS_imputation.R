library(dplyr)
library(ggplot2)

# Cargamos datos GSS 2004
gss_egor <- readRDS("trabajo_1_files/gss_egor.rds")
gss_egos <- gss_egor$ego
gss_alters <- gss_egor$alter

hist(gss_egos$age)
hist(gss_egos$educ_num)
unique(gss_egos$race)
unique(gss_egos$relig)

# Cargamos datos ATP W3

ATP_W3_sub <- readRDS("trabajo_1_files/ATP_W3_sub.rds")
labels <- sapply(ATP_W3_sub, function(x) attr(x, "label"))

library(haven)
ATP_W3 <- read_sav("B - Surveys Data/Datos A. Trends Panel/American-Trends-Panel-Wave-3-May-5-May-27/W3_May14/ATP W3.sav")

# Qué datos tenemos ---------------------------------

# Edad
labels[["F_AGECAT_TYPOLOGY"]]
ATP_W3$F_AGECAT_TYPOLOGY

# Educación
labels[["F_EDUCCAT_TYPOLOGY"]]
ATP_W3$F_EDUCCAT_TYPOLOGY

# Race
labels[["F_RACETHN_TYPOLOGY"]]
ATP_W3$F_RACETHN_TYPOLOGY

# Sex
labels[["F_SEX_FINAL"]]
ATP_W3$F_SEX_FINAL

# Religion
labels[["F_RELIG_TYPOLOGY"]]
ATP_W3$F_RELIG_TYPOLOGY

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
ggsave("trabajo_1_plots/age_GSS_vs_ATP.pdf", width = 8, height = 5)

# Comparación Chi-Cuadrado
# Necesitamos las cuentas, no las proporciones
gss_counts <- table(gss_egos$age_category_gss)
atp_counts_raw <- table(ATP_W3_sub$F_AGECAT_TYPOLOGY_factor)

# tabla de contingencia para el test
contingency_table_age <- rbind(GSS = as.numeric(gss_counts), ATP = as.numeric(atp_counts_raw))
colnames(contingency_table_age) <- age_labels_gss

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
set.seed(999)

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

pdf(file = "trabajo_1_plots/age_distribution_GSS.pdf", width = 6, height = 5)
hist(gss_egos$age, main = "Edades en GSS (EGO)", xlab = "Edad", ylab = "Frecuencia")
dev.off()

pdf(file = "trabajo_1_plots/age_distribution_ATP.pdf", width = 6, height = 5)
hist(ATP_W3_sub$age, main = "Edades Imputadas en ATP", xlab = "Edad", ylab = "Frecuencia")
dev.off()

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

# plot -----------
gss_counts <- table(gss_egos$degree_atp3cat)
atp_counts <- table(ATP_W3_sub$F_EDUCCAT_TYPOLOGY_factor)

gss_educ_proportions <- gss_counts / sum(gss_counts)
atp_educ_proportions <- atp_counts / sum(atp_counts)

df_gss_educ <- data.frame(
  source = "GSS",
  category = names(gss_educ_proportions), # O usa colnames(contingency_table_educ)
  proportion = as.numeric(gss_educ_proportions)
)

df_atp_educ <- data.frame(
  source = "ATP",
  category = names(atp_educ_proportions), # O usa colnames(contingency_table_educ)
  proportion = as.numeric(atp_educ_proportions)
)

educ_gss_atp <- rbind(df_gss_educ, df_atp_educ)

rm(df_gss_educ, df_atp_educ)

# Para ordenar los factores
educ_levels <- c("HS graduate or less", "Some college", "College graduate+")
educ_gss_atp$category <- factor(educ_gss_atp$category, levels = educ_levels)

print(educ_gss_atp)

plot_educ <- ggplot(educ_gss_atp, aes(x = category, y = proportion, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparación de Distribuciones de Nivel Educativo (GSS vs ATP)",
       x = "Nivel Educativo",
       y = "Proporción",
       fill = "Fuente de Datos") + # Etiqueta para la leyenda
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), # Puedes ajustar angle y hjust
        plot.title = element_text(hjust = 0.5)) # Centrar el título
print(plot_educ)
ggsave("trabajo_1_plots/educ_GSS_vs_ATP.pdf", plot = plot_educ, width = 8, height = 5)

# Tabla de contingencia
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

hist(gss_egos$educ_num, main = "Años de educación en GSS", xlab = "Años de educación")
hist(ATP_W3_sub$educ_num, main = "Años de educación imputados en ATP", xlab = "Años de educación")

pdf(file = "trabajo_1_plots/educ_distribution_GSS.pdf", width = 6, height = 5)
hist(gss_egos$educ_num, main = "Años de educación en GSS (EGO)", xlab = "Edad", ylab = "Frecuencia")
dev.off()

pdf(file = "trabajo_1_plots/educ_distribution_ATP.pdf", width = 6, height = 5)
hist(ATP_W3_sub$educ_num, main = "Años de educación imputados en ATP", xlab = "Edad", ylab = "Frecuencia")
dev.off()

# ------------------------------------------------------------------------------
# Raza
# ------------------------------------------------------------------------------


# --- Veamos cómo se comparan los niveles de RAZA en GSS y en ATP ---------

# Los niveles objetivo de GSS son: "Asian", "Black", "Hispanic", "White", "Other"
# ATP_W3_sub$F_RACETHN_TYPOLOGY:
# 1=White non-Hispanic, 2=Black non-Hispanic, 3=Hispanic, 4=Other, 9=DK/Ref

library(forcats)
gss_egos <- gss_egos %>%
  mutate(
    race_4cat = fct_collapse(race, Other = c("Asian", "Other")),
    race_4cat = factor(race_4cat, levels = c("Black", "Hispanic", "White", "Other"))
  )
print(table(gss_egos$race_4cat, useNA = "ifany"))

ATP_W3_sub <- ATP_W3_sub %>%
  mutate(
    race_harmonized = case_when(
      F_RACETHN_TYPOLOGY == 1 ~ "White",
      F_RACETHN_TYPOLOGY == 2 ~ "Black",
      F_RACETHN_TYPOLOGY == 3 ~ "Hispanic",
      F_RACETHN_TYPOLOGY == 4 ~ "Other", 
      F_RACETHN_TYPOLOGY == 9 ~ NA_character_,
      TRUE ~ NA_character_
    ),
    # Establecer los niveles para que coincidan con gss_egos$race_4cat
    race_harmonized = factor(race_harmonized, levels = c("Black", "Hispanic", "White", "Other"))
  )

# plot -----------
gss_race_counts_4cat <- table(gss_egos$race_4cat)
atp_race_counts_4cat <- table(ATP_W3_sub$race_harmonized)

gss_race_counts_4cat_proportions <- gss_race_counts_4cat / sum(gss_race_counts_4cat)
atp_race_counts_4cat_proportions <- atp_race_counts_4cat / sum(atp_race_counts_4cat)

df_gss_race <- data.frame(
  source = "GSS",
  category = names(gss_race_counts_4cat_proportions), # O usa colnames(contingency_table_educ)
  proportion = as.numeric(gss_race_counts_4cat_proportions)
)

df_atp_race <- data.frame(
  source = "ATP",
  category = names(atp_race_counts_4cat_proportions), # O usa colnames(contingency_table_educ)
  proportion = as.numeric(atp_race_counts_4cat_proportions)
)

race_gss_atp <- rbind(df_gss_race, df_atp_race)

rm(df_gss_race, df_atp_race)

# Para ordenar los factores
race_levels <- c("Black", "Hispanic", "White", "Other")
race_gss_atp$category <- factor(race_gss_atp$category, levels = race_levels)

print(race_gss_atp)

plot_educ <- ggplot(race_gss_atp, aes(x = category, y = proportion, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparación de Distribuciones de Raza (GSS vs ATP)",
       x = "Raza",
       y = "Proporción",
       fill = "Fuente de Datos") + # Etiqueta para la leyenda
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), # Puedes ajustar angle y hjust
        plot.title = element_text(hjust = 0.5)) # Centrar el título
print(plot_educ)
ggsave("trabajo_1_plots/race_GSS_vs_ATP.pdf", plot = plot_educ, width = 8, height = 5)

# Tabla de contingencia
contingency_table_race_4cat <- rbind(
  GSS = gss_race_counts_4cat[levels(gss_egos$race_4cat)],
  ATP = atp_race_counts_4cat[levels(ATP_W3_sub$race_harmonized)]
)
print(contingency_table_race_4cat)

# Test Chi-cuadrado
chi_sq_test_race_4cat <- chisq.test(contingency_table_race_4cat)

if (chi_sq_test_race_4cat$p.value < 0.05) {
  cat("Conclusión RAZA (4 cat): Las distribuciones de raza entre GSS y ATP son significativamente DIFERENTES (p < 0.05).\n")
} else {
  cat("Conclusión RAZA (4 cat): No hay evidencia de una diferencia significativa en las distribuciones de raza (p >= 0.05).\n")
}

plot(gss_egos$race_4cat, main = "Raza armonizada en GSS", xlab = "Raza")
plot(ATP_W3_sub$race_harmonized, main = "Raza armonizada en ATP", xlab = "Raza")


# ---------- Imputación de RAZA en ATP -----------------------------------------

gss_asian_other_subset <- gss_egos %>%
  filter(race %in% c("Asian", "Other")) # Usamos la variable original de 5 categorías

# proporciones de "Asian" y "Other" en este subconjunto
prob_asian <- prop.table(table(gss_asian_other_subset$race))["Asian"]
prob_other <- prop.table(table(gss_asian_other_subset$race))["Other"]

if (prob_asian + prob_other != 1) {cat("Cuidado, las probabilidades suman", prob_asian + prob_other, " ,y no 1")}

# nuevas columnas de raza imputada en ATP
ATP_W3_sub <- ATP_W3_sub %>%
  mutate(
    race_imputed_5cat = as.character(race_harmonized) # Empezar con las 4 categorías como character
  )
atp_other_indices <- which(ATP_W3_sub$race_harmonized == "Other") # índices de ATP que son "Other"

set.seed(456)

# Imputamos "Asian" u "Other" a los que tienen solo "Other" en ATP
for (i in atp_other_indices) {
  ATP_W3_sub$race_imputed_5cat[i] <- sample(
    x = c("Asian", "Other"),       
    size = 1,                      # Asignar una sola categoría
    replace = TRUE,
    prob = c(prob_asian, prob_other) # Con las probabilidades de GSS
  )
}
# Convertimos la nueva columna a factor con los 5 niveles
ATP_W3_sub$race_imputed_5cat <- factor(ATP_W3_sub$race_imputed_5cat,
                                       levels = c("Asian", "Black", "Hispanic", "White", "Other"))


# ------------------------------------------------------------------------------
# Raza
# ------------------------------------------------------------------------------

# --- Veamos cómo se comparan los niveles de RELIGIÓN en GSS y en ATP ----------

# Los niveles objetivo de GSS para religión son: "Catholic", "Jewish", "None", "OtherRelig", "Protestant"

ATP_W3_sub <- ATP_W3_sub %>%
  mutate(
    relig_harmonized = case_when(
      F_RELIG_TYPOLOGY == 1 ~ "Protestant", # ATP: 1 = Protestant (Baptist, Methodist, etc.)
      F_RELIG_TYPOLOGY == 2 ~ "Catholic",   # ATP: 2 = Roman Catholic (Catholic)
      F_RELIG_TYPOLOGY == 5 ~ "Jewish",     # ATP: 5 = Jewish (Judaism)
      F_RELIG_TYPOLOGY %in% c(9, 10, 12) ~ "None", # ATP: 9 = Atheist, 10 = Agnostic, 12 = Nothing in particular
      F_RELIG_TYPOLOGY %in% c(3, 4, 6, 7, 8, 11, 13, 14) ~ "OtherRelig",
                                            # ATP: 3 = Mormon
                                            # ATP: 4 = Orthodox (Greek, Russian, etc.)
                                            # ATP: 6 = Muslim (Islam)
                                            # ATP: 7 = Buddhist
                                            # ATP: 8 = Hindu
                                            # ATP: 11 = Something else (SPECIFY)
                                            # ATP: 13 = (VOL) Christian (no especificado como Protestante o Católico)
                                            # ATP: 14 = (VOL) Unitarian (Universalist)
      F_RELIG_TYPOLOGY == 99 ~ NA_character_, # ATP: 99 = (VOL) Don't know/Refused
      TRUE ~ NA_character_
    ),
    # Asegurar que los niveles del factor coincidan con los de gss_egos$relig
    relig_harmonized = factor(relig_harmonized, levels = c("Catholic", "Jewish", "None", "OtherRelig", "Protestant"))
  )

# Tabla de contingencia
gss_relig_counts <- table(gss_egos$relig)
atp_relig_counts <- table(ATP_W3_sub$relig_harmonized)
contingency_table_relig <- rbind(
  GSS = gss_relig_counts[levels(gss_egos$relig)],
  ATP = atp_relig_counts[levels(ATP_W3_sub$relig_harmonized)]
)
print(contingency_table_relig)

# Test Chi-cuadrado
chi_sq_test_relig <- chisq.test(contingency_table_relig)

if (chi_sq_test_relig$p.value < 0.05) {
  cat("Conclusión RELIGIÓN: Las distribuciones de religión entre GSS y ATP son significativamente DIFERENTES (p < 0.05).\n")
} else {
  cat("Conclusión RELIGIÓN: No hay evidencia de una diferencia significativa en las distribuciones de religión (p >= 0.05).\n")
}
