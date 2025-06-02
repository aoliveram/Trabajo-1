# -----------------------------------------------------------------------------
# ATP W3
# -----------------------------------------------------------------------------

library(haven)

ATP_W3 <- read_sav("B - Surveys Data/Datos A. Trends Panel/American-Trends-Panel-Wave-3-May-5-May-27/W3_May14/ATP W3.sav")
ATP_W3_df <- as.data.frame(ATP_W3)

# Extract labels
labels <- sapply(ATP_W3_df, function(x) attr(x, "label"))
labels[is.na(labels)] <- ""

# Convert the entire dataframe to character type
ATP_W3_df <- data.frame(lapply(ATP_W3_df, as.character), stringsAsFactors = FALSE)

# Add the labels as the first row
ATP_W3_df_2 <- rbind(labels, ATP_W3_df)

# Optionally, reset row names to avoid confusion
rownames(ATP_W3_df_2) <- NULL

# Write
write.csv(ATP_W3_df_2, file = "B - Surveys Data/Datos A. Trends Panel/ATP_W3.csv", row.names = FALSE)

# Exploración Variables -------------------------------------------------------

ATP_W3_df <- as.data.frame(ATP_W3)
labels <- sapply(ATP_W3_df , function(x) attr(x, "label"))

# Edad
labels[["F_AGECAT_TYPOLOGY"]]
ATP_W3_df$F_AGECAT_TYPOLOGY

# Educación
labels[["F_EDUCCAT_TYPOLOGY"]]
ATP_W3_df$F_EDUCCAT_TYPOLOGY

# Race
labels[["F_RACETHN_TYPOLOGY"]]
ATP_W3_df$F_RACETHN_TYPOLOGY

# Sex
labels[["F_SEX_FINAL"]]
ATP_W3_df$F_SEX_FINAL

# Religion
labels[["F_RELIG_TYPOLOGY"]]
ATP_W3_df$F_RELIG_TYPOLOGY

# Ver las primeras filas del data frame
head(ATP_W3)
head(ATP_W3_df)

# Calcular frecuencias y graficar las columnas seleccionadas
library(ggplot2)
library(gridExtra)

# Lista de columnas a graficar
columnas <- c("METECH_A_W3", "METECH_B_W3", "METECH_C_W3", 
              "METECH_D_W3", "METECH_E_W3", "METECH_F_W3")

# Calcular frecuencias para cada columna
frecuencias <- lapply(columnas, function(col) {
  table(factor(ATP_W3_df[[col]], levels = c(0, 1, 99)))
})

# Encontrar el valor máximo entre todas las frecuencias
max_y <- max(sapply(frecuencias, max))

# Crear los gráficos individuales
graficos <- lapply(1:length(columnas), function(i) {
  col <- columnas[i]
  frec <- as.data.frame(frecuencias[[i]])
  colnames(frec) <- c("Respuesta", "Frecuencia")
  
  ggplot(frec, aes(x = Respuesta, y = Frecuencia, fill = Respuesta)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("blue", "orange", "green")) +
    labs(title = col, x = "Respuestas", y = "Frecuencia") +
    ylim(0, max_y) +
    theme_minimal()
})

# Organizar los gráficos en un layout de 2 filas y 3 columnas
grid.arrange(grobs = graficos, nrow = 2)

############

delta_q = 1/6

# Descripciones de cada columna. 'early - st'
descripciones <- c(
  1 * delta_q, #"Usually try new products before others do",
  6 * delta_q, #"Prefer my tried and trusted brands",
  3 * delta_q, #"Like being able to tell others about new brands and products I have tried",
  2 * delta_q, #"Like the variety of trying new products",
  5 * delta_q, #"Feel more comfortable using familiar brands and products",
  4 * delta_q #"Wait until I hear about others’ experiences before I try new products"
)

# Crear un data frame para almacenar los resultados
resultados <- data.frame(
  Columna = columnas,
  Descripcion = descripciones,
  Total_Respuestas_0_1 = numeric(length(columnas)),
  Porcentaje_Respuestas_1 = numeric(length(columnas))
)

# Calcular los valores para cada columna
for (i in seq_along(columnas)) {
  col <- columnas[i]
  
  # Frecuencias de la columna
  frec <- table(factor(ATP_W3_df[[col]], levels = c(0, 1, 99)))
  
  # Sumar respuestas 0 y 1
  total_0_1 <- sum(frec[c("0", "1")])
  
  # Porcentaje de respuestas 1 con respecto al total de respuestas (0 + 1)
  porcentaje_1 <- (frec["1"] / total_0_1) * 100
  
  # Guardar resultados en el data frame
  resultados$Total_Respuestas_0_1[i] <- total_0_1
  resultados$Porcentaje_Respuestas_1[i] <- porcentaje_1
}

# Mostrar los resultados
print(resultados)

# Calcular la suma total de respuestas 1 en todas las columnas
suma_respuestas_1 <- sum(ATP_W3_df[columnas] == 1, na.rm = TRUE)

# Mostrar el resultado
print(suma_respuestas_1)


# ----------------------------- Cramos .rds -----------------------------------

library(haven) 
library(dplyr)

ATP_W3 <- read_sav("B - Surveys Data/Datos A. Trends Panel/American-Trends-Panel-Wave-3-May-5-May-27/W3_May14/ATP W3.sav")
ATP_W3_df <- as.data.frame(ATP_W3)

# --- Columnas Demográficas ATP ---
columnas_demograficas_atp <- c(
  "QKEY",                  # Identificador único del encuestado
  "F_AGECAT_TYPOLOGY",     # Edad
  "F_EDUCCAT_TYPOLOGY",    # Educación
  "F_RACETHN_TYPOLOGY",    # Raza
  "F_SEX_FINAL",           # Sexo
  "F_RELIG_TYPOLOGY",      # Religión
  "WEIGHT_W3"              # La variable de peso
)

# --- Columnas de Innovación ATP ---
columnas_innovacion_atp_w3 <- c(
  "METECH_A_W3", "METECH_B_W3", "METECH_C_W3",
  "METECH_D_W3", "METECH_E_W3", "METECH_F_W3",
  "MEFOOD_A_W3", "MEFOOD_B_W3", "MEFOOD_C_W3",
  "MEFOOD_D_W3", "MEFOOD_E_W3", "MEFOOD_F_W3"
)

# --- Combinamos ---
columnas_seleccionadas_total_w3 <- c(columnas_demograficas_atp, columnas_innovacion_atp_w3)

# --- Creamos DataFrame ---
ATP_W3_sub <- ATP_W3_df %>%
  select(all_of(columnas_seleccionadas_total_w3))

# --- Guardar el DataFrame de subconjunto ---

saveRDS(ATP_W3_sub, file = "trabajo_1_files/ATP_W3_sub.rds")

# Convertir columnas que son factores (o etiquetadas por haven) a character para CSV
ATP_W3_sub_csv <- ATP_W3_sub
for (col_name in names(ATP_W3_sub_csv)) {
  if (is.factor(ATP_W3_sub_csv[[col_name]]) || inherits(ATP_W3_sub_csv[[col_name]], "haven_labelled")) {
    ATP_W3_sub_csv[[col_name]] <- as.character(haven::as_factor(ATP_W3[[col_name]])) # Usar as_factor para obtener etiquetas
  }
}
# También convierte las variables METECH a character si aún son etiquetadas
for (col_name in columnas_innovacion_atp) {
  if (col_name %in% names(ATP_W3_sub_csv) && inherits(ATP_W3[[col_name]], "haven_labelled")){
    ATP_W3_sub_csv[[col_name]] <- as.character(haven::as_factor(ATP_W3[[col_name]]))
  }
}

write.csv(ATP_W3_sub_csv, file = "trabajo_1_files/ATP_W3_sub.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# ATP W4
# -----------------------------------------------------------------------------

library(haven)

ATP_W4 <- read_sav("B - Surveys Data/Datos A. Trends Panel/American-Trends-Panel-Wave-4-May-30-Jun-30/W4_Jun14/ATP W4.sav") # Leer el archivo como tibble
ATP_W4_df <- as.data.frame(ATP_W4)

# Extract labels
labels <- sapply(ATP_W4_df, function(x) attr(x, "label"))
labels[is.na(labels)] <- ""

# Convert the entire dataframe to character type
ATP_W4_df <- data.frame(lapply(ATP_W4_df, as.character), stringsAsFactors = FALSE)

# Add the labels as the first row
ATP_W4_df_2 <- rbind(labels, ATP_W4_df)

# Optionally, reset row names to avoid confusion
rownames(ATP_W4_df_2) <- NULL

# Write
write.csv(ATP_W4_df_2, file = "B - Surveys Data/Datos A. Trends Panel/ATP_W4.csv", row.names = FALSE)

# Exploración Variables -------------------------------------------------------

ATP_W4_df <- as.data.frame(ATP_W4)
labels <- sapply(ATP_W4_df , function(x) attr(x, "label"))

# Edad
labels[["F_AGECAT_TYPOLOGY"]]
ATP_W4_df$F_AGECAT_TYPOLOGY

hist(as.numeric(ATP_W4_df$F_AGECAT_TYPOLOGY),
     breaks = 0:4,
     xaxt = "n",
     main = "Distribución por Categoría de Edad",
     xlab = "Rango de Edad")
axis(1, at = 1:4, labels = c("18-29", "30-49", "50-64", "65+"))


# Gráfico de barras
barplot(table(factor(
                ATP_W4_df$F_AGECAT_TYPOLOGY,
                levels = 1:4,
                labels = c("18-29", "30-49", "50-64", "65+")
              )),
        col = "skyblue",
        xlab = "Rango de Edad",
        ylab = "Frecuencia",
        main = "Distribución por Categoría de Edad")


# Educación
labels[["F_EDUCCAT_TYPOLOGY"]]
ATP_W4_df$F_EDUCCAT_TYPOLOGY

# Race
labels[["F_RACETHN_TYPOLOGY"]]
ATP_W4_df$F_RACETHN_TYPOLOGY

# Sex
labels[["F_SEX_FINAL"]]
ATP_W4_df$F_SEX_FINAL

# Religion
labels[["F_RELIG_TYPOLOGY"]]
ATP_W4_df$F_RELIG_TYPOLOGY

# Ver las primeras filas del data frame
head(ATP_W4)
head(ATP_W4_df)

# ----------------------------- Cramos .rds -----------------------------------

library(haven) 
library(dplyr)

# --- Columnas Demográficas ATP ---
columnas_demograficas_atp <- c(
  "QKEY",                  # Identificador único del encuestado
  "F_AGECAT_TYPOLOGY",     # Edad
  "F_EDUCCAT_TYPOLOGY",    # Educación
  "F_RACETHN_TYPOLOGY",    # Raza
  "F_SEX_FINAL",           # Sexo
  "F_RELIG_TYPOLOGY",      # Religión
  "WEIGHT_W4"              # La variable de peso
)

# --- Columnas de Innovación ATP ---
columnas_innovacion_atp_w4 <- c(
  "METECH_A_W4", "METECH_B_W4", "METECH_C_W4",
  "METECH_D_W4", "METECH_E_W4", "METECH_F_W4",
  "MEFOOD_A_W4", "MEFOOD_B_W4", "MEFOOD_C_W4",
  "MEFOOD_D_W4", "MEFOOD_E_W4", "MEFOOD_F_W4"
)

# --- Combinamos ---
columnas_seleccionadas_total_w4 <- c(columnas_demograficas_atp, columnas_innovacion_atp_w4)

# --- Creamos DataFrame ---
ATP_W4_sub <- ATP_W4_df %>%
  select(all_of(columnas_seleccionadas_total_w4))

# --- Guardar el DataFrame de subconjunto ---

saveRDS(ATP_W4_sub, file = "trabajo_1_files/ATP_W4_sub.rds")

# Convertir columnas que son factores (o etiquetadas por haven) a character para CSV
ATP_W4_sub_csv <- ATP_W4_sub
for (col_name in names(columnas_demograficas_atp)) {
  if (is.factor(ATP_W4_sub_csv[[col_name]]) || inherits(ATP_W4_sub_csv[[col_name]], "haven_labelled")) {
    ATP_W4_sub_csv[[col_name]] <- as.character(haven::as_factor(ATP_W4[[col_name]])) # Usar as_factor para obtener etiquetas
  }
}
# También convierte las variables METECH a character si aún son etiquetadas
for (col_name in columnas_innovacion_atp) {
  if (col_name %in% names(ATP_W4_sub_csv) && inherits(ATP_W4[[col_name]], "haven_labelled")){
    ATP_W4_sub_csv[[col_name]] <- as.character(haven::as_factor(ATP_W4[[col_name]]))
  }
}

write.csv(ATP_W4_sub_csv, file = "trabajo_1_files/ATP_W4_sub.csv", row.names = FALSE)
