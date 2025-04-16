# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Objetivo: Crear un objeto EGOR a partir de datos GSS 2004 y simular redes
# Datos: GSS_2004_EGO_bin.csv
# Paquetes: egor, ergm.ego, dplyr, tidyr, readr
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- 1. Carga de Paquetes ---
# Instalar paquetes si no están presentes
# install.packages(c("egor", "ergm.ego", "dplyr", "tidyr", "readr", "srvyr"))
library(egor)
library(ergm.ego)
library(dplyr)
library(tidyr)
library(readr)
library(srvyr) # Necesario para el objeto egor aunque no se use directamente aquí

# --- 2. Carga de Datos ---
gss_data_raw <- read_csv("trabajo_1_files/GSS_2004_EGO_bin.csv", na = "NA", show_col_types = FALSE)

# Inspección inicial
# glimpse(gss_data_raw)
# head(gss_data_raw)
# summary(gss_data_raw$numgiven)

# --- 3. Preparación de Datos para EGOR ---

# 3.1. Crear DataFrame de Egos (egos)
# Añadir un ID único para cada ego
egos_df <- gss_data_raw %>%
  mutate(.egoID = row_number()) %>% # Crear ID único para cada ego
  select(
    .egoID, # ID único del ego
    sex, race, educ, age, relig, degree, # Atributos del ego
    numgiven # Número de alters nominados
  )

# Revisar tipos de datos y recodificar si es necesario (ejemplo para factores)
# egos_df <- egos_df %>%
#   mutate(
#     sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
#     race = factor(race, levels = c(1, 2, 3), labels = c("White", "Black", "Other"))
#     # Añadir más recodificaciones según sea necesario basándose en el codebook
#   )

# glimpse(egos_df)
# head(egos_df)

# 3.2. Crear DataFrame de Alters (alters)
# Necesitamos transformar los datos de formato ancho a largo

# Columnas de atributos de los alters
alter_attr_cols <- c(
  paste0("sex", 1:5), paste0("race", 1:5), paste0("educ", 1:5),
  paste0("age", 1:5), paste0("relig", 1:5)
)

# Columnas de estado/relación de los alters
alter_status_cols <- c(
  paste0("spouse", 1:5), paste0("parent", 1:5), paste0("sibling", 1:5),
  paste0("child", 1:5), paste0("othfam", 1:5), paste0("cowork", 1:5),
  paste0("memgrp", 1:5), paste0("neighbr", 1:5), paste0("friend", 1:5),
  paste0("advisor", 1:5), paste0("other", 1:5), paste0("talkto", 1:5)
)

# Combinar todas las columnas de alters
alter_cols_all <- c(alter_attr_cols, alter_status_cols)

# Crear el DataFrame de alters en formato largo
alters_long <- gss_data_raw %>%
  mutate(.egoID = row_number()) %>% # Añadir ID de ego para unir/filtrar
  select(.egoID, numgiven, all_of(alter_cols_all)) %>%
  pivot_longer(
    cols = all_of(alter_cols_all),
    names_to = c(".value", ".alterID_char"), # .value captura el nombre base, .alterID_char el número
    names_pattern = "([a-zA-Z]+)(\\d+)", # Separa nombre (letras) de número (dígitos)
    values_drop_na = FALSE # Mantener filas aunque el valor original sea NA (importante antes de filtrar por numgiven)
  ) %>%
  mutate(.alterID = as.integer(.alterID_char)) %>% # Convertir ID de alter a número
  filter(.alterID <= numgiven) %>% # Mantener solo alters nominados
  select(-numgiven, -.alterID_char) # Eliminar columnas auxiliares

# Revisar la estructura y datos de alters
# glimpse(alters_long)
# head(alters_long)
# summary(alters_long)

# Nota sobre codificación:
# Los atributos como 'educ1', 'relig1', 'race1' tienen códigos específicos del GSS.
# Para usarlos en modelos (ej. nodematch), podrían necesitar recodificación a categorías más simples
# o ser tratados como numéricos/ordinales si tiene sentido.
# Las variables de estado (spouse1, parent1, etc.) son 1=Yes, 2=No. Para usarlas como indicadores
# binarios, podríamos recodificarlas a 1=Yes, 0=No. O crear una única variable categórica 'relationship'.
# Por ahora, las mantenemos como están, pero esto es importante para el modelado.

# Ejemplo de recodificación para variables de estado (opcional, depende del modelo):
# alters_long <- alters_long %>%
#   mutate(across(starts_with(c("spouse", "parent", "sibling", "child", "othfam",
#                              "cowork", "memgrp", "neighbr", "friend", "advisor", "other")),
#                 ~ if_else(. == 1, 1, 0, missing = 0))) # 1 si es Yes, 0 si No o NA

# 3.3. Crear DataFrame de Lazos Alter-Alter (aaties) - CORREGIDO
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

# Revisar la estructura y datos de aaties (opcional)
# glimpse(aaties_long)
# nrow(aaties_long) # Verificar cuántos lazos se encontraron

head(alters_long)
head(egos_df)
head(aaties_long)

# Crear el objeto egor
# Puede que necesites ajustar los nombres en ID.vars si los cambiaste arriba
my_gss_egor <- egor(
  alters = alters_long,
  egos = egos_df,
  aaties = aaties_long,
  ID.vars = list(
    ego = ".egoID",    # Columna de ID del ego en todos los dataframes
    alter = ".alterID", # Columna de ID del alter (dentro del ego) en alters_long
    source = ".srcID", # Columna de ID del alter origen en aaties_long
    target = ".tgtID"  # Columna de ID del alter destino en aaties_long
  ),
  alter_design = list(max = 5)
)

# Verificar el objeto creado
print(my_gss_egor)
summary(my_gss_egor)

# --- 5. Preparación para Modelado con ergm.ego ---

# 5.1. Recodificación de Variables (Ejemplos - ¡Ajustar según necesidad!)
# Es *crucial* recodificar las variables categóricas/binarias adecuadamente antes de modelar.

# Recodificar sexo (Ego y Alter) a factor
my_gss_egor$ego <- my_gss_egor$ego %>%
  mutate(sex_cat = factor(sex, levels = c(1, 2), labels = c("Male", "Female")))
my_gss_egor$alter <- my_gss_egor$alter %>%
  mutate(sex_cat = factor(sex, levels = c(1, 2), labels = c("Male", "Female")))

# Recodificar raza (Ego y Alter) a factor (agrupando 'Other')
my_gss_egor$ego <- my_gss_egor$ego %>%
  mutate(race_cat = factor(race, levels = c(1, 2, 3), labels = c("White", "Black", "Other")))
my_gss_egor$alter <- my_gss_egor$alter %>%
  mutate(race_cat = case_when(
    race == 1 ~ "Asian", # Basado en codebook RACE1
    race == 2 ~ "Black",
    race == 3 ~ "Hispanic",
    race == 4 ~ "White",
    race == 5 ~ "Other",
    TRUE ~ NA_character_ # Maneja NA, 9, etc.
  )) %>%
  mutate(race_cat = factor(race_cat)) # Convertir a factor

# Recodificar variables de estado (relación) a binario 0/1
status_vars <- c("spouse", "parent", "sibling", "child", "othfam",
                 "cowork", "memgrp", "neighbr", "friend", "advisor", "other")
my_gss_egor$alter <- my_gss_egor$alter %>%
  mutate(across(all_of(status_vars), ~ if_else(. == 1, 1, 0, missing = 0), .names = "{.col}_bin"))

# Recodificar educ (Alter) basado en codebook EDUC1 (ejemplo agrupando)
# 0: <HS, 1: HS, 2: Some College/Assoc, 3: Bach, 4: Grad
my_gss_egor$alter <- my_gss_egor$alter %>%
  mutate(educ_cat = case_when(
    educ %in% c(0, 1, 2) ~ "Less than HS", # Codes 0, 1, 2
    educ == 3           ~ "High School",  # Code 3
    educ %in% c(4, 5)   ~ "Some College", # Codes 4, 5
    educ == 6           ~ "Bachelor",     # Code 6
    educ == 7           ~ "Graduate",     # Code 7
    TRUE                ~ NA_character_   # Codes 8, 9 (DK/NA) or others
  )) %>%
  mutate(educ_cat = factor(educ_cat, levels = c("Less than HS", "High School", "Some College", "Bachelor", "Graduate")))

# Tratar edad (Alter) como numérica (ya debería serlo, pero verificamos)
my_gss_egor$alter <- my_gss_egor$alter %>% mutate(age = as.numeric(age))

# Verificar las nuevas variables
# glimpse(my_gss_egor$ego)
# glimpse(my_gss_egor$alter)
# summary(my_gss_egor$alter$race_cat)
# summary(my_gss_egor$alter$educ_cat)

# --- 6. Definición y Ajuste del Modelo ERGM-EGO ---

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


# 6.1. Definir la fórmula del modelo SIMPLE
# Modelo MUY simple: Solo tendencia a formar lazos y efecto del sexo del ego

# Verificar si hay lazos alter-alter (necesario para 'edges')
n_aaties <- nrow(my_gss_egor$aatie)
cat("\nNúmero de lazos alter-alter encontrados:", n_aaties, "\n")

# 6.2. Ajustar el modelo (si hay una fórmula válida)

# Usar un tamaño de muestra pequeño para probar/rapidez inicial
sample_size_fit <- min(nrow(my_gss_egor$ego), 200) # Reducido para prueba rápida
# Para el análisis final, usar: sample_size_fit <- nrow(my_gss_egor$ego)

cat("\nAjustando el modelo ERGM-EGO SIMPLIFICADO con la fórmula:\n")
print(formula_to_use)
cat("Usando tamaño de muestra:", sample_size_fit, "\n")

# ¡Esto debería ser mucho más rápido!

fitted_model <- ergm.ego(
  my_gss_egor ~ edges
)

head(my_gss_egor$ego)
head(my_gss_egor$alter)
unique(my_gss_egor$ego$race_cat)
head(my_gss_egor$aatie)

fitted_model <- ergm.ego(
  my_gss_egor ~ edges + degree(0), 
  #control = snctrl(ppopsize=1000)
)

# --- 7. Simulación de Redes ---
simulated_networks <- NULL
if (!is.null(fitted_model)) {
  n_sim <- 1000 # Número de redes a simular
  cat("\nSimulando", n_sim, "redes desde el modelo ajustado...\n")
  
  # La simulación también puede ser intensiva
  simulated_networks <- simulate(
    fitted_model,
    nsim = n_sim
    # constraints = ~ blockdiag(".egoID") # Asegura que los lazos solo ocurren dentro de cada ego-red
    # Podrían necesitarse otros parámetros/controles dependiendo del modelo
  )
  
  cat("Simulación completada.\n")
  
  # Inspeccionar las redes simuladas (opcional)
  # print(simulated_networks) # Es una lista de objetos network/igraph
  # plot(simulated_networks[[1]]) # Graficar la primera red simulada
  # summary(simulated_networks[[1]] ~ edges + nodematch("sex_cat")) # Calcular estadísticas en una red simulada
  
} else {
  cat("\nNo se pueden simular redes porque el modelo no fue ajustado.\n")
}

# --- Fin del Script ---
cat("\nScript finalizado.\n")
# El objeto 'my_gss_egor' está listo para análisis.
# Si el modelo se ajustó, 'fitted_model' contiene los resultados.
# Si la simulación se realizó, 'simulated_networks' contiene la lista de redes simuladas.