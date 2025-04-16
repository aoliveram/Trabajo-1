# --- 5. Preparación para Modelado con ergm.ego (CORREGIDO Y AMPLIADO) ---

# Recodificar sexo (Ego y Alter) a factor (ya hecho, verificar consistencia)
my_gss_egor$ego <- my_gss_egor$ego %>%
  mutate(sex_cat = factor(sex, levels = c(1, 2), labels = c("Male", "Female")))
my_gss_egor$alter <- my_gss_egor$alter %>%
  mutate(sex_cat = factor(sex, levels = c(1, 2), labels = c("Male", "Female")))

# --- ARMONIZAR Raza ---
# Mapear alter race a las categorías del ego
my_gss_egor$alter <- my_gss_egor$alter %>%
  mutate(race_cat_harmonized = case_when(
    race == 4 ~ "White",    # Code 4 = White
    race == 2 ~ "Black",    # Code 2 = Black
    race %in% c(1, 3, 5) ~ "Other", # Codes 1(Asian), 3(Hisp), 5(Other) -> Other
    TRUE ~ NA_character_
  )) %>%
  # Asegurar que los niveles del factor sean los mismos que en ego$race_cat
  mutate(race_cat_harmonized = factor(race_cat_harmonized, levels = c("White", "Black", "Other")))

# Recodificar ego race (ya hecho, verificar consistencia)
my_gss_egor$ego <- my_gss_egor$ego %>%
  mutate(race_cat_harmonized = factor(race, levels = c(1, 2, 3), labels = c("White", "Black", "Other")))

# Recodificar educ (Alter) basado en codebook EDUC1 (como antes)
my_gss_egor$alter <- my_gss_egor$alter %>%
  mutate(educ_cat = case_when(
    educ %in% c(0, 1, 2) ~ "Less than HS",
    educ == 3           ~ "High School",
    educ %in% c(4, 5)   ~ "Some College",
    educ == 6           ~ "Bachelor",
    educ == 7           ~ "Graduate",
    TRUE                ~ NA_character_
  )) %>%
  mutate(educ_cat = factor(educ_cat, levels = c("Less than HS", "High School", "Some College", "Bachelor", "Graduate")))

# Tratar edad (Alter) como numérica (ya hecho)
my_gss_egor$alter <- my_gss_egor$alter %>% mutate(age = as.numeric(age))

# --- Convertir Variables de Estado Binarias a Factor ---
status_vars_to_factor <- c("spouse", "parent", "sibling", "child", "othfam",
                           "cowork", "memgrp", "neighbr", "friend", "advisor", "other")

my_gss_egor$alter <- my_gss_egor$alter %>%
  mutate(across(all_of(status_vars_to_factor),
                ~ factor(if_else(. == 1, "Yes", "No", missing = "No")), # Crear factores Yes/No
                .names = "{.col}_fct")) # Nuevos nombres para no sobreescribir _bin si se crearon

# Verificar las nuevas variables y niveles
# glimpse(my_gss_egor$ego)
# glimpse(my_gss_egor$alter)
# levels(my_gss_egor$ego$race_cat)
# levels(my_gss_egor$alter$race_cat_harmonized) # Deberían coincidir
# class(my_gss_egor$alter$age)
# class(my_gss_egor$alter$spouse_fct)
# class(my_gss_egor$alter$friend_fct)


# --- 6. Definición y Ajuste del Modelo ERGM-EGO (Revisado) ---

# 6.1. Definir la fórmula del modelo REVISADA
# Usar las variables armonizadas y los términos correctos (nodecov/nodefactor)

n_aaties <- nrow(my_gss_egor$aatie)
cat("\nNúmero de lazos alter-alter encontrados:", n_aaties, "\n")

formula_to_use <- NULL

if (n_aaties > 50) { # Mantener umbral (ajustar si es necesario)
  # Fórmula revisada
  model_formula_net_revised <- ~ edges +
    nodematch("sex_cat", diff = TRUE) +
    nodematch("race_cat_harmonized", diff = TRUE) + # Usar la variable armonizada
    nodematch("educ_cat", diff = TRUE) +
    nodecov("age") +                  # Usar nodecov para edad numérica
    nodefactor("spouse_fct") +        # Usar nodefactor para el factor spouse
    nodefactor("friend_fct") +        # Usar nodefactor para el factor friend
    ego.nodefactor("sex_cat") +       # Esto usa ego$sex_cat
    ego.nodefactor("race_cat")        # Añadir efecto de la raza del ego
  
  formula_to_use <- model_formula_net_revised
  cat("Usando fórmula revisada con términos de red alter-alter:\n")
  print(formula_to_use)
  
} else if (n_aaties > 0) {
  # Mantener el modelo simple si hay pocos lazos
  simple_formula <- ~ edges + ego.nodefactor("sex_cat")
  formula_to_use <- simple_formula
  cat("Pocos lazos alter-alter. Usando fórmula simple solo con 'edges'.\n")
  print(formula_to_use)
  
}


head(my_gss_egor$ego)
head(my_gss_egor$alter)
unique(my_gss_egor$ego$race_cat)
head(my_gss_egor$aatie)

fitted_model <- ergm.ego(
  my_gss_egor ~ edges + nodematch("sex_cat", diff = TRUE)
  #control = snctrl(ppopsize=1000)
)
