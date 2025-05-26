library(network)
library(dplyr)
library(ggplot2)
#library(haven) # Para el manejo inicial si los datos vienen de .sav
library(gridExtra)


# --- 1. Cargamos datos ---

load("trabajo_1_files/ATP_network_simulated_1000.RData")

# --- 2. Extraer los atributos METECH_ de la red a un dataframe para facilitar el cálculo ---
metech_vars <- c("metech_a", "metech_b", "metech_c", 
                 "metech_d", "metech_e", "metech_f")

df_metech_attr <- data.frame(QKEY = 1:network.size(ATP_network_simulated_1000)) # Crear un ID temporal

for (m_var in metech_vars) {
  if (m_var %in% list.vertex.attributes(ATP_network_simulated_1000)) {
    df_metech_attr[[m_var]] <- get.vertex.attribute(ATP_network_simulated_1000, m_var)
  } else {
    warning(paste("Atributo", m_var, "no encontrado en la red. Se creará con NAs."))
    df_metech_attr[[m_var]] <- NA
  }
}

# --- 3. Distribución original METECH ---

plots_metech_original <- list()
for (i in seq_along(metech_vars)) {
  col_name <- metech_vars[i]
  
  # Manejar NAs para la tabla y el gráfico
  valid_data <- na.omit(df_metech_attr[[col_name]])
  if (length(valid_data) == 0) {
    cat("No hay datos válidos para", col_name, "\n")
    next
  }
  
  freq_table <- as.data.frame(table(factor(valid_data, levels = c(0, 1))))
  colnames(freq_table) <- c("Respuesta", "Frecuencia")
  
  p <- ggplot(freq_table, aes(x = Respuesta, y = Frecuencia, fill = Respuesta)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("0" = "skyblue", "1" = "coral")) + # Colores para 0 y 1
    labs(title = paste("Var:", col_name), x = "(0=No, 1=Sí)", y = "Frecuencia") +
    ylim(-1,1000) +
    theme_minimal() +
    theme(legend.position = "none")
  plots_metech_original[[col_name]] <- p
}

pdf("trabajo_1_plots/metech_distribution_ATP.pdf", width = 10, height = 7)
do.call(grid.arrange, c(plots_metech_original, ncol = 3))
dev.off()
