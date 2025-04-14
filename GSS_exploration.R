# ----------------- DESCRIPTORES ----------------- 
# 2812 obs , 1264 var --> El 50.7% de ellos tienen entrada no nula en 'numgiven'.
# 1426 obs con 'numgiven'
# 1070 obs con 'numgiven' != 0

# ----------------- IMPORTAR ------------

library(haven)

GSS_2004 <- read_dta("B - Surveys Data/GSS 2004/GSS 2004 NORC.dta")

class(GSS_2004)
# Entradas únicas: 0, 1, 2, 3, 4, 5, 6. No hay '9'.
unique(GSS_2004$numgiven)
class(GSS_2004$numgiven[1])
# Entradas enteras
sum(GSS_2004$numgiven %% 1 == 0, na.rm = TRUE) / nrow(GSS_2004)
# Son equivalentes: NA == NA(i)
sum(is_tagged_na(GSS_2004$numgiven))
sum(is.na(GSS_2004$numgiven))
sum(!is.na(GSS_2004$numgiven))

# Veamos los lazos alter-alter
unique(GSS_2004$close12)
plot(GSS_2004$close12)

# ----------------- PROCESAMIENTO ------------

library(haven)
library(dplyr)
library(labelled)

# 1. Importar datos
GSS_2004 <- read_dta("B - Surveys Data/GSS 2004/GSS 2004 NORC.dta")

# 2. Definir columnas de atributos y de red
col_alters_attr <- c(
  "numgiven", 
  paste0("sex", 1:5), paste0("race", 1:5), paste0("educ", 1:5), 
  paste0("age", 1:5), paste0("relig", 1:5)
)
col_alters_net <- paste0("close", c("12", "13", "14", "15", "23", "24", "25", "34", "35", "45"))

# 3. Filtrar solo casos con numgiven no nulo
GSS_2004_EGO <- GSS_2004 %>% filter(!is.na(numgiven))

# 4. Seleccionar columnas relevantes (atributos + red)
selected_columns <- c(col_alters_attr, col_alters_net)
GSS_2004_EGO <- GSS_2004_EGO %>% select(all_of(selected_columns))

# 5. Guardar etiquetas originales
labels_list_attr <- lapply(GSS_2004_EGO[col_alters_attr], val_labels)
labels_list_net  <- lapply(GSS_2004_EGO[col_alters_net],  val_labels)

# 6. Convertir a data frame base R 
GSS_2004_EGO <- GSS_2004_EGO %>% mutate(across(everything(), ~ as.integer(.)))
GSS_2004_EGO <- as.data.frame(GSS_2004_EGO)

# 7. Análisis exploratorio

library(ggplot2)

# Resumen estadístico 
summary(GSS_2004_EGO)

# Valores faltantes por columna
missing_counts <- colSums(is.na(GSS_2004_EGO))
print(missing_counts)
plot(missing_counts)

# Distribución de numgiven
ggplot(GSS_2004_EGO, aes(x = numgiven)) +
  geom_bar(fill = "darkorange", color = "black", alpha = 0.7) +
  labs(title = "Número de personas mencionadas (numgiven)", x = "numgiven", y = "Frecuencia") +
  theme_minimal()

sum(GSS_2004_EGO$numgiven != 0, na.rm = TRUE) # 1070 personas reportaron cercanos

# Histogramas de todas las variables (excepto numgiven)
plot_histograms <- function(df) {
  num_vars <- setdiff(names(df), col_alters_net)  # Excluir numgiven si necesario
  
  plots <- lapply(num_vars, function(var) {
    ggplot(df, aes(x = .data[[var]])) +
      geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.7) +
      labs(title = paste("Distribución de", var), x = var, y = "Frecuencia") +
      theme_minimal()
  })
  
  return(plots)
}
plot_histograms(GSS_2004_EGO)

# 8. Exportamos datos limpios
write.csv(GSS_2004_EGO, 'B - Surveys Data/GSS 2004/GSS_2004_EGO.csv', row.names = FALSE)
GSS_2004_EGO <- read.csv('B - Surveys Data/GSS 2004/GSS_2004_EGO.csv')

# ----------------- NETWORK ------------

