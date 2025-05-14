library(dplyr)

gss_egor <- readRDS("trabajo_1_files/gss_egor.rds")
gss_egos <- gss_egor$ego
gss_alters <- gss_egor$alter

# --- Crear educ_level_unified para Egos ---
# Los niveles de gss_egos$degree son: LTHS, HS, Assoc, Bach, Grad
# Los niveles de gss_alters$educ_cat son: LTHS, HS, SomeColl, Assoc, Bach, Grad
# Unificaremos a: "LTHS", "HS", "SomeCollege/Assoc", "Bach", "Grad"

gss_egos <- gss_egos %>%
  mutate(
    educ_level_unified = factor(case_when(
      degree == "LTHS" ~ "LTHS",
      degree == "HS" ~ "HS",
      degree == "Assoc" ~ "SomeCollege/Assoc", # Mapea Assoc a la categoría combinada
      degree == "Bach" ~ "Bach",
      degree == "Grad" ~ "Grad",
      TRUE ~ NA_character_ # Por si hay algún NA o valor inesperado en degree
    ), levels = c("LTHS", "HS", "SomeCollege/Assoc", "Bach", "Grad"))
  )

# --- Crear educ_level_unified para Alters ---
# Aquí, combinamos "SomeColl" y "Assoc" en una sola categoría.
# La variable 'educ' en gss_alters es la que usaste para crear educ_cat.
# educ_cat = case_when(
#       educ %in% c(0, 1, 2) ~ "LTHS",
#       educ == 3 ~ "HS",
#       educ == 4 ~ "SomeColl",
#       educ == 5 ~ "Assoc",
#       educ == 6 ~ "Bach",
#       educ == 7 ~ "Grad",
#       TRUE ~ NA_character_
# )
# Así que gss_alters$educ_cat ya existe y tiene los niveles deseados,
# solo necesitamos re-factorizar para combinar "SomeColl" y "Assoc".

gss_alters <- gss_alters %>%
  mutate(
    educ_level_unified_temp = educ_cat, # Usamos la columna existente
    educ_level_unified = factor(
      ifelse(educ_level_unified_temp %in% c("SomeColl", "Assoc"), 
             "SomeCollege/Assoc", 
             as.character(educ_level_unified_temp)
      ),
      levels = c("LTHS", "HS", "SomeCollege/Assoc", "Bach", "Grad")
    )
  ) %>% select(-educ_level_unified_temp) # Quitamos la columna temporal

# Verificar los niveles
print("Niveles de educ_level_unified en egos:")
print(table(gss_egos$educ_level_unified, useNA = "ifany"))
print("Niveles de educ_level_unified en alters:")
print(table(gss_alters$educ_level_unified, useNA = "ifany"))

# Ahora, actualiza el objeto egor con estas nuevas columnas
# Es importante que los dataframes dentro de gss_egor sean los actualizados.
gss_egor$ego <- gss_egos
gss_egor$alter <- gss_alters
