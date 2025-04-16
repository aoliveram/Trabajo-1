# --- 0. Load Necessary Libraries ---
# install.packages(c("readr", "dplyr", "tidyr", "egor", "ergm.ego", "network"))
library(readr)
library(dplyr)
library(tidyr)
library(egor)
library(ergm.ego)
library(network)

# --- 1. Load the Raw Data ---
gss_data_raw <- read_csv("trabajo_1_files/GSS_2004_EGO_bin.csv", na = "NA", show_col_types = FALSE)

# Add a unique Ego ID
gss_data_raw <- gss_data_raw %>%
  mutate(.egoID = row_number()) %>%
  relocate(.egoID)

# --- 2. Prepare Ego Data ------------------------------------------------------
# Keep ALL egos, including those with numgiven = 0
egos_df <- gss_data_raw %>%
  select(
    .egoID,
    sex, race, educ, age, relig, degree, numgiven # Keep numgiven temporarily for alter processing
  ) %>%
  # --- Recode Ego variables ---
  mutate(
    # Ensure NAs are handled *before* factoring
    sex = na_if(sex, 9), # Assuming 9 might be an NA code not mentioned
    race = na_if(race, 9), # Assuming 9 might be an NA code
    relig = na_if(relig, 98), # Codebook NA/DK
    relig = na_if(relig, 99),
    degree = na_if(degree, 8), # Codebook DK/NA
    degree = na_if(degree, 9),
    numgiven = na_if(numgiven, 9), # Codebook NA for numgiven is 9
    
    # Sex: 1=Male, 2=Female
    sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    
    # Race: Harmonize to White, Black, OtherRace
    race = case_when(
      race == 1 ~ "White",
      race == 2 ~ "Black",
      race == 3 ~ "OtherRace", # Ego code 3 is 'Other'
      TRUE ~ NA_character_
    ),
    race = factor(race, levels = c("White", "Black", "OtherRace")),
    
    # Educ: Numeric years. Code 99 = NA.
    educ_num = ifelse(educ == 99, NA_real_, as.numeric(educ)),
    
    # Age: Numeric. Assuming direct read is fine, add specific NA checks if needed.
    age = as.numeric(age),
    
    # Relig: Simplify categories
    relig = case_when(
      relig == 1 ~ "Protestant",
      relig == 2 ~ "Catholic",
      relig == 3 ~ "Jewish",
      relig == 4 ~ "None",
      relig %in% c(5, 6, 7, 8, 9, 10, 11, 12, 13) ~ "OtherRelig",
      TRUE ~ NA_character_ # Catches NAs assigned above and any others
    ),
    relig = factor(relig),
    
    # Degree: Simplify categories
    degree = case_when(
      degree == 0 ~ "LTHS",
      degree == 1 ~ "HS",
      degree == 2 ~ "Assoc",
      degree == 3 ~ "Bach",
      degree == 4 ~ "Grad",
      TRUE ~ NA_character_ # Catches NAs assigned above
    ),
    degree = factor(degree, levels = c("LTHS", "HS", "Assoc", "Bach", "Grad"))
  ) %>% 
  select(-educ) # Remove original numeric code

# --- 3. Prepare Alter Data ----------------------------------------------------

alter_attr_vars <- c("sex", "race", "educ", "age", "relig")
alter_status_vars <- c("spouse", "parent", "sibling", "child", "othfam",
                       "cowork", "memgrp", "neighbr", "friend", "advisor",
                       "other", "talkto")

alter_attr_cols <- paste0(alter_attr_vars, rep(1:5, each = length(alter_attr_vars)))
alter_status_cols <- paste0(alter_status_vars, rep(1:5, each = length(alter_status_vars)))

# Pivot longer first
alters_long <- gss_data_raw %>%
  select(.egoID, numgiven, all_of(alter_attr_cols), all_of(alter_status_cols)) %>%
  # Handle potential NA in numgiven *before* filtering
  mutate(numgiven = na_if(numgiven, 9)) %>%
  pivot_longer(
    cols = c(all_of(alter_attr_cols), all_of(alter_status_cols)),
    names_to = c(".variable", "alt_ID"),
    names_pattern = "([a-zA-Z]+)(\\d+)",
    values_to = ".value",
    values_drop_na = FALSE # Keep all rows for now
  ) %>%
  mutate(alt_ID = as.integer(alt_ID))

# Filter *only* based on whether the alter position was nominated
# Keep egos with numgiven=0 (they will have no rows here after the filter)
# Keep egos with numgiven=NA (they will also have no rows)
alters_long_filtered <- alters_long %>%
  filter(!is.na(numgiven) & alt_ID <= numgiven) # Keep only valid alter slots

# Pivot wider
alters_wide <- alters_long_filtered %>%
  pivot_wider(
    names_from = .variable,
    values_from = .value
  ) %>%
  select(-numgiven) # numgiven no longer needed in alters table

# --- Recode Alter variables ---
alters_df <- alters_wide %>%
  mutate(
    # Ensure NAs are handled *before* factoring for Alters
    sex = na_if(sex, 9),
    race = na_if(race, 9),
    educ = na_if(educ, 8), # DK
    educ = na_if(educ, 9), # NA
    relig = na_if(relig, 8), # DK
    relig = na_if(relig, 9), # NA
    # Handle NAs for all status vars (using 9)
    across(all_of(alter_status_vars), ~na_if(., 9)),
    
    # Sex: 1=Male, 2=Female
    sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    
    # Race: Harmonize to White, Black, OtherRace (matching ego structure)
    race = case_when(
      race == 4 ~ "White",     # Alter code 4 = White
      race == 2 ~ "Black",     # Alter code 2 = Black
      race %in% c(1, 3, 5) ~ "OtherRace", # Alter codes 1=Asian, 3=Hispanic, 5=Other -> OtherRace
      TRUE ~ NA_character_ # Catches NAs assigned above and any others
    ),
    # IMPORTANT: Ensure levels match ego$race levels exactly
    race = factor(race, levels = c("White", "Black", "OtherRace")),
    
    # Educ: Map to simpler categories & create numeric version
    # Alter codes: 0=1-6, 1=7-9, 2=10-12, 3=HSG, 4=SomeColl, 5=Assoc, 6=Bach, 7=Grad
    educ_cat = case_when(
      educ %in% c(0, 1, 2) ~ "LTHS",
      educ == 3 ~ "HS",
      educ == 4 ~ "SomeColl",
      educ == 5 ~ "Assoc",
      educ == 6 ~ "Bach",
      educ == 7 ~ "Grad",
      TRUE ~ NA_character_ # Catches NAs assigned above
    ),
    educ_cat = factor(educ_cat, levels = c("LTHS", "HS", "SomeColl", "Assoc", "Bach", "Grad")),
    # Numeric version - Approximation (can refine midpoints if needed)
    educ_num = case_when(
      educ == 0 ~ 4, # Midpoint approx for 1-6
      educ == 1 ~ 8, # Midpoint approx for 7-9
      educ == 2 ~ 11, # Midpoint approx for 10-12
      educ == 3 ~ 12, # HS Grad
      educ == 4 ~ 14, # Some College approx
      educ == 5 ~ 14, # Assoc approx
      educ == 6 ~ 16, # Bach
      educ == 7 ~ 18, # Grad approx
      TRUE ~ NA_real_
    ),
    
    # Age: Numeric
    age = as.numeric(age),
    
    # Relig: Simplify categories
    relig = case_when(
      relig == 1 ~ "Protestant",
      relig == 2 ~ "Catholic",
      relig == 3 ~ "Jewish",
      relig %in% c(4, 8) ~ "None",
      relig == 5 ~ "OtherRelig",
      TRUE ~ NA_character_ # Catches NAs assigned above
    ),
    relig = factor(relig),
    
    # Relationship Statuses: 1=Yes, 2=No -> Logical
    across(
      all_of(alter_status_vars),
      ~ case_when(
        . == 1 ~ TRUE,
        . == 2 ~ FALSE,
        TRUE ~ NA      # Includes NAs assigned above and any others
      )
    )
  ) %>%
  # Filter on NAs in sex and race
  filter(!is.na(sex), !is.na(race)) %>%
  # Keep categorical educ, rename numeric for clarity if using both
  select(-educ) # Remove original numeric code

# Remove numgiven from egos_df now that alter processing is done
egos_df <- egos_df %>% select(-numgiven)

# --- Display structure and summary ---

cols_common <- c('.egoID', 'sex', 'race', 'age', 'relig', 'educ_num')

# Verifica que ambas tengan las columnas requeridas
cat("¿egos_df tiene todas las columnas? ", all(cols_common %in% names(egos_df)), "\n")
cat("¿alters_df tiene todas las columnas? ", all(cols_common %in% names(alters_df)), "\n\n")

# Compara valores únicos por columna
lapply(cols_common, function(col) {
  same_vals <- setequal(unique(egos_df[[col]]), unique(alters_df[[col]]))
  class_ego <- class(egos_df[[col]])
  class_alter <- class(alters_df[[col]])
  
  cat("Columna:", col, "\n")
  cat("  ¿Valores únicos iguales?:", same_vals, "\n")
  cat("  Clase en egos_df:", class_ego, "\n")
  cat("  Clase en alters_df:", class_alter, "\n\n")
})

# --- 3.3. DataFrame Lazos Alter-Alter (aaties) --------------------------------
# Columnas de la red alter-alter
alter_net_cols <- paste0("close", c("12", "13", "14", "15", "23", "24", "25", "34", "35", "45"))

# Crear el DataFrame de lazos en formato largo (edge list)
aaties_long <- gss_data_raw %>%
  mutate(.egoID = row_number()) %>%
  # Seleccionar .egoID, numgiven y las columnas de lazos alter-alter
  select(.egoID, numgiven, all_of(alter_net_cols)) %>%
  pivot_longer(
    cols = all_of(alter_net_cols),
    names_to = "connection",
    values_to = "tie_value",
    values_drop_na = FALSE # Es importante no descartar NA antes de filtrar por tie_value == 1
  ) %>%
  # Filtrar solo las conexiones existentes (donde closeXY == 1)
  filter(tie_value == 1) %>%
  # Extraer los IDs de origen y destino del nombre de la conexión (closeXY)
  mutate(
    .srcID = as.integer(substr(connection, nchar("close") + 1, nchar("close") + 1)),
    .tgtID = as.integer(substr(connection, nchar("close") + 2, nchar("close") + 2))
  ) %>%
  # Asegurarse de que ambos alters conectados estén dentro del numgiven del ego
  # La columna 'numgiven' ya está presente desde el select inicial
  filter(.srcID <= numgiven & .tgtID <= numgiven) %>%
  # Seleccionar solo las columnas necesarias para el formato de aaties
  select(.egoID, .srcID, .tgtID)


# --- 4. Create the EGOR Object ------------------------------------------------

head(egos_df)
head(alters_df)
head(aaties_long)

sum(is.na(egos_df$sex))
sum(is.na(egos_df$race))
sum(is.na(egos_df$age))
sum(is.na(egos_df$relig))
sum(is.na(egos_df$educ_num))

# alters_df <- alters_df[!is.na(egos_df$age), ] # quitamos primero en alters
# egos_df <- egos_df[!is.na(egos_df$age), ]
# alters_df <- alters_df[!is.na(egos_df$relig), ]
# egos_df <- egos_df[!is.na(egos_df$relig), ]
# alters_df <- alters_df[!is.na(egos_df$educ_num), ]
# egos_df <- egos_df[!is.na(egos_df$educ_num), ]

sum(is.na(alters_df$sex))
sum(is.na(alters_df$race))
sum(is.na(alters_df$age))
sum(is.na(alters_df$relig))
sum(is.na(alters_df$educ_num))

# alters_df <- alters_df[!is.na(alters_df$age), ]
# alters_df <- alters_df[!is.na(alters_df$relig), ]
# alters_df <- alters_df[!is.na(alters_df$educ_num), ]

# 1. Filtra egos
egos_df <- egos_df[complete.cases(egos_df[, c("sex", "race", "age", "relig", "educ_num")]), ]

# 2. Filtra alters por egos válidos
alters_df <- alters_df[alters_df$.egoID %in% egos_df$.egoID, ]

# 3. (Opcional) Filtra alters por NA en variables de interés
alters_df <- alters_df[complete.cases(alters_df[, c("sex", "race", "age", "relig", "educ_num")]), ]

# 4. Filtra aaties por egos y alters válidos
aaties_long <- aaties_long[
  aaties_long$.egoID %in% egos_df$.egoID &
    aaties_long$.srcID %in% alters_df$alt_ID &
    aaties_long$.tgtID %in% alters_df$alt_ID, ]

# 1. Crea un vector con los alters válidos para cada ego
valid_alters <- alters_df[, c(".egoID", "alt_ID")]

# 2. Une (merge) para verificar que tanto .srcID como .tgtID existen para el ego correspondiente
aaties_long <- merge(
  aaties_long, valid_alters, 
  by.x = c(".egoID", ".srcID"), by.y = c(".egoID", "alt_ID")
)
aaties_long <- merge(
  aaties_long, valid_alters, 
  by.x = c(".egoID", ".tgtID"), by.y = c(".egoID", "alt_ID")
)

# 3. (Opcional) Si quieres dejar solo las columnas originales:
aaties_long <- aaties_long[, c(".egoID", ".srcID", ".tgtID")]



# Access data within the object using $ego and $alter
gss_egor <- egor(
  egos = egos_df,
  alters = alters_df,
  aaties = aaties_long,
  ID.vars = list(
    ego = ".egoID",    # Columna de ID del ego en todos los dataframes
    alter = "alt_ID", # Columna de ID del alter (dentro del ego) en alters_long
    source = ".srcID", # Columna de ID del alter origen en aaties_long
    target = ".tgtID"  # Columna de ID del alter destino en aaties_long
  ),
  alter_design = list(max = 5)
)

summary(gss_egor)

head(gss_egor$ego)
head(gss_egor$alter)
head(gss_egor$aatie)

# --- 5. Define and Estimate ERGM.EGO Model ---

# - Tendencia general a formar lazos alter-alter (edges)
# - Homofilia por sexo entre alters (nodematch('sex_cat'))
# - Homofilia por raza entre alters (nodematch('race_cat'))
# - Efecto del sexo del ego en el número de lazos alter-alter (ego.nodefactor('sex_cat'))
# - Efecto de la edad del alter en la probabilidad de tener lazos (nodefactor('age')) - ¡Necesitaría discretizar o usar nodecov!
# - Tendencia de cónyuges a estar conectados con otros alters (nodefactor('spouse_bin'))

# NOTA IMPORTANTE: La presencia de muy pocos lazos alter-alter (aaties) puede hacer
# que los modelos con términos de red alter-alter (como edges, nodematch) sean
# difíciles o imposibles de estimar (degeneración). Si aaties_long está vacío o casi vacío,
# estos términos no deben incluirse.


# `nodefactor` requiere que no hayan NA.

ergm.ego(gss_egor ~ edges + nodematch("race") + nodematch("sex") )
         # control = control.ergm.ego(MPLE.samplesize = 10000), # Example control if needed

ergm.ego(gss_egor ~ edges + nodematch("race") + nodecov("age"))




#
# --- Analizar el modelo y perfeccionarlo ---------------
# ~
# ~
# ~
# ~
# -------------------------------------------------------





# --- 6. Check Fit and Simulate Networks ---
if (!is.null(fit_gss_ego)) {
  cat("\n--- Model Estimation Summary ---\n")
  print(summary(fit_gss_ego))
  
  # Check MCMC diagnostics if needed (requires coda package)
  # mcmc.diagnostics(fit_gss_ego)
  
  cat("\nSimulating 1000 networks based on the fitted model...\n")
  simulated_networks <- simulate(fit_gss_ego, nsim = 1000, verbose = TRUE)
  
  cat(sprintf("\nSimulation complete. %d sets of ego networks were generated.\n", length(simulated_networks)))
  
  # Example: Get summary statistics from simulated networks
  # This calculates the average degree for each ego across all simulations
  # avg_sim_degrees <- sapply(seq_len(nrow(gss_egor$ego)), function(i) {
  #    mean(sapply(simulated_networks, function(net_list) network.edgecount(net_list[[i]])))
  # })
  # print("Average simulated degrees for first 10 egos:")
  # print(head(avg_sim_degrees, 10))
  
  # Save results
  # saveRDS(fit_gss_ego, file = "gss_ego_model_fit.rds")
  # saveRDS(simulated_networks, file = "gss_simulated_networks.rds")
  
} else {
  cat("\nERGM.EGO model estimation failed. Cannot proceed with simulation.\n")
  cat("Potential issues:\n")
  cat("- Model misspecification (terms might be collinear or poorly identified).\n")
  cat("- Insufficient MCMC convergence (try adjusting control.ergm.ego parameters).\n")
  cat("- Remaining data inconsistencies (double-check factor levels and NAs).\n")
}

cat("\nScript finished.\n")