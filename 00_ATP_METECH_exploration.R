library(haven)

# Leer el archivo .sav y convertirlo en un data frame
#ruta_archivo <- "C:/Users/Usuario/Desktop/Datos A. Trends Panel/American-Trends-Panel-Wave-3-May-5-May-27/W3_May14/ATP W3.sav"#"ruta/del/archivo.sav"  # Cambia esto por la ruta real de tu archivo

ATP_W3 <- read_sav("B - Surveys Data and Papers Inn/Datos A. Trends Panel/American-Trends-Panel-Wave-3-May-5-May-27/W3_May14/ATP W3.sav") # Leer el archivo como tibble
ATP_W3_df <- as.data.frame(ATP_W3)        # Convertir a data frame (opcional si prefieres este formato)

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
