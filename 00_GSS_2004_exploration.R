# ----------------- DESCRIPTORES ----------------- 
# 2812 obs , 1264 var --> El 50.7% de ellos tienen entrada no nula en 'numgiven'.
# 1426 obs con 'numgiven' != NA
# 1070 obs con 'numgiven' != 0
# 789 obs con 'numgiven' >= 2

# NO tomé una decisión con close_recode: 
# Está NA --> 0 (no lazo), pero NO hay NA en GSS_2004_EGO

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

# 1) Importar datos
GSS_2004 <- read_dta("B - Surveys Data/GSS 2004/GSS 2004 NORC.dta")

# 1.1) Variable de Pesos

GSS_2004$wtss
GSS_2004$wtssall
GSS_2004$wtssnr

max(unique(GSS_2004$wtssnr))
min(unique(GSS_2004$wtssnr))

weights_2004 <- GSS_2004$wtssnr

1413 * min(unique(weights_2004)) / 3

# 2) Definir columnas de atributos y de red

col_ego_attr <- c("sex", "race", "educ", "age", "relig",  "degree", "wtssnr")

col_alters_attr <- c(
  "numgiven", 
  paste0("sex", 1:5), paste0("race", 1:5), paste0("educ", 1:5), 
  paste0("age", 1:5), paste0("relig", 1:5)
)

col_alters_net <- paste0("close", c("12", "13", "14", "15", "23", "24", "25", "34", "35", "45"))

col_alters_status <- c(
  paste0("spouse", 1:5),
  paste0("parent", 1:5),
  paste0("sibling", 1:5),
  paste0("child", 1:5),
  paste0("othfam", 1:5),
  paste0("cowork", 1:5),
  paste0("memgrp", 1:5),
  paste0("neighbr", 1:5),
  paste0("friend", 1:5),
  paste0("advisor", 1:5),
  paste0("other", 1:5),
  paste0("talkto", 1:5)
)

# 3) Filtrar solo casos con numgiven no nulo
GSS_2004_EGO <- GSS_2004 %>% filter(!is.na(numgiven))

# 4) Seleccionar columnas relevantes (atributos + red)
selected_columns <- c(col_ego_attr, col_alters_attr, col_alters_net, col_alters_status)
GSS_2004_EGO <- GSS_2004_EGO %>% select(all_of(selected_columns))

# 5) Guardar etiquetas originales
labels_ego_attr <- lapply(GSS_2004_EGO[col_alters_attr], val_labels)
labels_alters_attr <- lapply(GSS_2004_EGO[col_alters_attr], val_labels)
labels_alters_net  <- lapply(GSS_2004_EGO[col_alters_net],  val_labels)
labels_alters_status <- lapply(GSS_2004_EGO[col_alters_status], val_labels)

# 6) Convertir a data frame base R 
#GSS_2004_EGO <- GSS_2004_EGO %>% mutate(across(everything(), ~ as.integer(.)))
GSS_2004_EGO <- GSS_2004_EGO %>%
  mutate(
    across(-c(wtssnr), ~ as.integer(.)),
    wtssnr = as.numeric(wtssnr)
  )
GSS_2004_EGO <- as.data.frame(GSS_2004_EGO)

# 7) Análisis exploratorio

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

# 8) Exportamos datos limpios
write.csv(GSS_2004_EGO, 'trabajo_1_files/GSS_2004_EGO.csv', row.names = FALSE)
GSS_2004_EGO <- read.csv('trabajo_1_files/GSS_2004_EGO.csv')

saveRDS(GSS_2004_EGO, 'trabajo_1_files/GSS_2004_EGO.rds')

# ----------------- NETWORK ------------

# 9) Cuantificar casos '7 = refused'
refused_count <- sum(sapply(GSS_2004_EGO[col_alters_net], function(x) sum(x == 7, na.rm = TRUE)))
cat("Cantidad de casos '7 = refused':", refused_count, "\n")

# 10) Recodificar lazos alter-alter y guardar
set.seed(123) # Para reproducibilidad

recode_close <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) return(0)           # NA --> 0 (no hay NA en GSS_2004_EGO)
    if (val == 3) return(0)             # 'total strangers' --> 0
    if (val %in% c(1, 2)) return(1)     # 'especially close' 'know each other' --> 1
    if (val == 7) return(sample(c(0, 1), 1)) # 'refused' --> random (0,1)
    return(NA)
  })
}

GSS_2004_EGO_bin <- GSS_2004_EGO
GSS_2004_EGO_bin[col_alters_net] <- lapply(GSS_2004_EGO[col_alters_net], recode_close)

write.csv(GSS_2004_EGO_bin, 'trabajo_1_files/GSS_2004_EGO_bin.csv', row.names = FALSE)
GSS_2004_EGO <- read.csv('trabajo_1_files/GSS_2004_EGO_bin.csv')

saveRDS(GSS_2004_EGO_bin, 'trabajo_1_files/GSS_2004_EGO_bin.rds')

# 11) Contar el número de egos para cada valor de alters declarados
table_numgiven <- table(GSS_2004_EGO_bin$numgiven)

percent_numgiven <- table_numgiven / sum(table_numgiven)

png("trabajo_1_plots/n_alters_plot.png", width = 800, height = 600)
n_alters_plot <- barplot(percent_numgiven,
              col = "skyblue",
              border = "black",
              main = "Egos que declararon n Alters - GSS_2004_EGO_bin",
              xlab = "Número de Alters (n)",
              ylab = "Porcentaje de Egos",
              names.arg = 0:6,
              ylim = c(0, max(percent_numgiven) * 1.1))
text(x = n_alters_plot, 
     y = percent_numgiven, 
     labels = as.numeric(table_numgiven), 
     pos = 3, cex = 1, col = "black")
text(x = max(n_alters_plot), 
     y = max(percent_numgiven) * 0.85, 
     labels = paste("Total de casos:", sum(table_numgiven)), 
     adj = 1, cex = 1.1, font = 2)
dev.off()

# 12) Graficar la red alter-alter para un ego

library(igraph)

plot_ego_network <- function(row, close_cols) {
  adj <- matrix(0, nrow = 5, ncol = 5)
  colnames(adj) <- rownames(adj) <- paste0("Alter", 1:5)
  alter_pairs <- list(
    c(1,2), c(1,3), c(1,4), c(1,5),
    c(2,3), c(2,4), c(2,5),
    c(3,4), c(3,5),
    c(4,5)
  )
  for (i in seq_along(close_cols)) {
    val <- as.numeric(row[[close_cols[i]]])
    if (!is.na(val) && val == 1) {
      a <- alter_pairs[[i]][1]
      b <- alter_pairs[[i]][2]
      adj[a, b] <- 1
      adj[b, a] <- 1
    }
  }
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  plot(g, main = "Red alter-alter (un ego)", vertex.size = 30)
}

# Graficamos redes únicas
plot_ego_network(GSS_2004_EGO_bin[3, ], col_alters_net)

# 13) Histograma de densidad de redes Alter-Alter (lazos máximos 5x4/2 = 10)
densidades <- rowSums(GSS_2004_EGO_bin[, col_alters_net], na.rm = TRUE) / 10

png("trabajo_1_plots/densidades.png", width = 800, height = 600)
hist(densidades,
     breaks = seq(0, 1, by = 0.1),
     col = "cyan",
     border = "black",
     main = "Densidades de red (basado 5 alters)",
     xlab = "Densidad",
     ylab = "Frecuencia",
     xlim = c(0, 1))
dev.off()

# redes con densidad mayor a 0
sum(densidades > 0)

# 14) Histograma con densidad corregida por Alter reportados. Eje-y considera número de redes con al menos 2 Alters.

calc_density <- function(row, close_cols) {
  n <- row[["numgiven"]]
  if (n < 2) {#print('NA')
    return(NA)} # No se puede calcular densidad con menos de 2 alters 
                # Ojo que choose(0, 2)=choose(1, 2)=0, se indetermina la frac --> NA
  
  # Calculamos número máximo de links n(n-1)/2 .
  possible_links <- choose(n, 2)
  # Extrae los valores de los lazos relevantes. Ej: close12, close13, close23
  links <- as.numeric(row[close_cols][1:possible_links])
  observed_links <- sum(links, na.rm = TRUE)
  
  return(observed_links / possible_links)
}

# Aplica la función a cada fila
densidades_corregidas <- apply(GSS_2004_EGO_bin, 1, calc_density, close_cols = col_alters_net)

# Filtra densidades válidas (no NA)
densidades_validas <- densidades_corregidas[!is.na(densidades_corregidas)] 
length(densidades_validas)/length(densidades_corregidas)
length(densidades_validas) # 789 obs con 'numgiven' >= 2

# Histograma con eje-y ajustado
png("trabajo_1_plots/densidades_corregidas.png", width = 800, height = 600)
hist(densidades_validas,
     breaks = seq(0, 1, by = 0.1),
     col = "cyan",
     border = "black",
     main = "Densidades de red (corregido por redes con $n>=2$ )",
     xlab = "Densidad corregida",
     ylab = "Frecuencia",
     xlim = c(0, 1),
     ylim = c(0, length(densidades_validas)))
dev.off
