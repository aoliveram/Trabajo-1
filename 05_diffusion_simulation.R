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

# --- 4. Índice de Propensión a la Innovación ---
# a. Usually try new products before others do (METECH_A_W3): 1 = pro-innovación
# b. Prefer my tried and trusted brands (METECH_B_W3): 1 = ANTI-innovación (invertir)
# c. Like being able to tell others about new brands and products I have tried (METECH_C_W3): 1 = pro-innovación
# d. Like the variety of trying new products (METECH_D_W3): 1 = pro-innovación
# e. Feel more comfortable using familiar brands and products (METECH_E_W3): 1 = ANTI-innovación (invertir)
# f. Wait until I hear about others' experiences before I try new products (METECH_F_W3): 1 = ANTI-innovación (invertir)

df_metech_attr <- df_metech_attr %>%
  mutate(
    # Si el valor original es NA, el resultado de la recodificación también será NA
    score_A = ifelse(is.na(metech_a), NA, metech_a),                             # Pro
    score_B = ifelse(is.na(metech_b), NA, 1 - metech_b),                         # ANTI
    score_C = ifelse(is.na(metech_c), NA, metech_c),                             # Pro
    score_D = ifelse(is.na(metech_d), NA, metech_d),                             # Pro
    score_E = ifelse(is.na(metech_e), NA, 1 - metech_e),                         # ANTI
    score_F = ifelse(is.na(metech_f), NA, 1 - metech_f))                         # ANTI

df_metech_attr <- df_metech_attr %>%
  mutate(
    # sumamos los scores (0-6). # si hay entrada NA en alguna, el resultado es NA.
    new_tech_pref_raw_index = rowSums(select(., starts_with("score_"))),
    
    # Normalizar el índice para que esté entre 0 y 1
    alpha_innov_prop = new_tech_pref_raw_index / 6
  )

# --- 5. Estudiar Distribución del Índice Creado y Guardar Histograma ---
cat("\nResumen del índice alpha_innov_prop:\n")
print(summary(df_metech_attr$alpha_innov_prop))
cat("\nTabla de frecuencia para el índice crudo (0-6):\n")
print(table(df_metech_attr$new_tech_pref_raw_index, useNA = "ifany"))


# Histograma del índice normalizado (0-1)
hist_alpha <- ggplot(df_metech_attr, aes(x = alpha_innov_prop)) +
  geom_histogram(binwidth = 1/6 / 2, fill = "darkgreen", color = "black", boundary = 0) + # binwidth ajustado
  scale_x_continuous(breaks = seq(0, 1, by = 1/6), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Distribución del Índice de Propensión a la Innovación (alpha_innov_prop)",
       x = "Índice Normalizado (0-1)", y = "Frecuencia") +
  theme_minimal()
print(hist_alpha)
ggsave("trabajo_1_plots/metech_alpha_distribution_ATP.pdf", plot = hist_alpha, width = 8, height = 6)
