library(network)
library(dplyr)
library(ggplot2)
library(gridExtra)

# --- Variables y parámetros iniciales ---
metech_vars <- c("metech_a", "metech_b", "metech_c", 
                 "metech_d", "metech_e", "metech_f")

networks_dir <- "trabajo_1_files/ATP_network_ergm/"
output_prefix <- "ATP_network_simulated_1000_mur_"
N_networks <- 100 # Cambiar si tienes un número distinto de redes

# --- Loop sobre las redes simuladas ---
for (i in 1:N_networks) {
  
  # Cargar red
  filename <- sprintf("%sATP_network_simulated_1000_%03d.rds", networks_dir, i)
  ATP_network_simulated_1000 <- readRDS(filename)
  
  # --- 2. Extraer atributos METECH ---
  df_metech_attr <- data.frame(QKEY = 1:network.size(ATP_network_simulated_1000)) # ID temporal
  
  for (m_var in metech_vars) {
    if (m_var %in% list.vertex.attributes(ATP_network_simulated_1000)) {
      df_metech_attr[[m_var]] <- get.vertex.attribute(ATP_network_simulated_1000, m_var)
    } else {
      warning(paste("Atributo", m_var, "no encontrado en la red", i, ". Se creará con NAs."))
      df_metech_attr[[m_var]] <- NA
    }
  }
  
  # --- 4. Índice de Propensión a la Innovación ---
  df_metech_attr <- df_metech_attr %>%
    mutate(
      score_A = ifelse(is.na(metech_a), NA, metech_a),
      score_B = ifelse(is.na(metech_b), NA, 1 - metech_b),
      score_C = ifelse(is.na(metech_c), NA, metech_c),
      score_D = ifelse(is.na(metech_d), NA, metech_d),
      score_E = ifelse(is.na(metech_e), NA, 1 - metech_e),
      score_F = ifelse(is.na(metech_f), NA, 1 - metech_f)
    ) %>%
    mutate(
      new_tech_pref_raw_index = rowSums(select(., starts_with("score_"))),
      q_i = new_tech_pref_raw_index / 6
    )
  
  # --- 5. Asignar atributo a la red ---
  set.vertex.attribute(ATP_network_simulated_1000, "q_i", df_metech_attr$q_i)
  
  # --- 6. Guardar nueva red con el atributo alpha ---
  output_filename <- sprintf("%s%s%03d.rds", networks_dir, output_prefix, i)
  saveRDS(ATP_network_simulated_1000, output_filename)
}


# -----------------------------------------------------------------------------
# Para solo la red ATP simulada original
# -----------------------------------------------------------------------------


# --- 1. Cargamos datos ---

ATP_network_simulated_1000 <- readRDS("trabajo_1_files/ATP_network_simulated_1000.rds")

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

# --- 5. Distribución del Índice ---
print(table(df_metech_attr$new_tech_pref_raw_index, useNA = "ifany"))

# Crear una columna factor para graficar, y asegurar que los NAs se manejen
df_metech_attr <- df_metech_attr %>% 
  mutate(alpha_innov_prop_factor = factor(alpha_innov_prop,
                                          levels = seq(0, 1, by = 1/6),
                                          labels = round(seq(0, 1, by = 1/6),3)))

# Gráfico de barras corregido usando geom_bar y el factor
hist_alpha <- ggplot(df_metech_attr, aes(x = alpha_innov_prop_factor)) +
  geom_bar(fill = "green", color = "black", stat = "count") + # stat="count" es el default para geom_bar
  labs(title = "Distribución del Índice de Propensión a la Innovación (alpha_innov_prop)",
       x = "Alpha Index", 
       y = "Frecuencia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) # Ajusta el ángulo si es necesario

print(hist_alpha)
ggsave("trabajo_1_plots/metech_alpha_distribution_ATP.pdf", plot = hist_alpha, width = 8, height = 6)
got
# --- 6. Agregarmos índice como atributo a los nodos de la red ---
set.vertex.attribute(ATP_network_simulated_1000, 
                     "alpha_innov_prop", 
                     df_metech_attr$alpha_innov_prop)

print(list.vertex.attributes(ATP_network_simulated_1000))

saveRDS(ATP_network_simulated_1000, "trabajo_1_files/ATP_network_simulated_1000_alpha.rds")
ATP_network_simulated_1000_alpha <- readRDS("trabajo_1_files/ATP_network_simulated_1000_alpha.rds")
