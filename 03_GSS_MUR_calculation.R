# ------------------------------------------------------------------------------
# Script: 05_mur_GSS_calculation.R
#
# Objetivo:
#   Calcula y asigna el 'Score de Propensión a la Acción Colectiva Normativa'
#   a cada nodo en las 100 redes GSS simuladas. El score se basa en 9 variables
#   de comportamiento y disposición política.
#
# Modifica:
#   - Procesa los archivos 'GSS_network_simulated_1000_XXX.rds'.
#   - Guarda nuevas versiones de las redes con el score añadido, usando el
#     prefijo '_mur_' (ej. 'GSS_network_simulated_1000_mur_XXX.rds').
#
# Genera (pdf):
#   - trabajo_1_plots/gss_propensity_vars_distribution.pdf
#   - trabajo_1_plots/gss_propensity_score_distribution.pdf
# ------------------------------------------------------------------------------

# --- 0. Cargar Librerías ---
library(network)
library(dplyr)
library(ggplot2)
library(gridExtra)

# --- 1. Variables y Parámetros Iniciales ---

# Las 9 variables "ingrediente" para nuestro score
propensity_ingredient_vars <- c("signdpet", "avoidbuy", "joindem", "attrally", 
                                "cntctgov", "polfunds", "usemedia", "interpol", "actlaw")

# Rutas y nombres de archivos
networks_dir <- "trabajo_1_files/GSS_network_ergm/"
output_prefix <- "GSS_network_simulated_1000_mur_"
N_networks <- 100

# --- 2. Loop sobre las 100 Redes Simuladas ---

cat(paste("Procesando", N_networks, "redes GSS para calcular el score de propensión...\n"))

for (i in 1:N_networks) {
  
  # Construir el nombre del archivo de entrada
  input_filename <- file.path(networks_dir, sprintf("GSS_network_simulated_1000_%03d.rds", i))
  
  if (!file.exists(input_filename)) {
    warning(paste("Archivo no encontrado, saltando:", input_filename))
    next
  }
  
  # Cargar la red simulada
  gss_network <- readRDS(input_filename)
  
  # Extraer los atributos de propensión a un dataframe para facilitar el cálculo
  df_propensity_attr <- data.frame(.ID = 1:network.size(gss_network)) # ID temporal
  
  for (p_var in propensity_ingredient_vars) {
    if (p_var %in% list.vertex.attributes(gss_network)) {
      df_propensity_attr[[p_var]] <- get.vertex.attribute(gss_network, p_var)
    } else {
      warning(paste("Atributo", p_var, "no encontrado en la red", i))
      df_propensity_attr[[p_var]] <- NA
    }
  }
  
  # --- 3. Calcular el Score de Propensión a la Acción Colectiva ---
  
  # Primero, recodificar las variables a una escala común (0 a 3) donde 3 es alta propensión
  df_propensity_attr <- df_propensity_attr %>%
    mutate(
      # Para las primeras 8 variables (1=hecho reciente -> 3, 4=nunca -> 0)
      recod_signdpet = case_when(signdpet == 1 ~ 3, signdpet == 2 ~ 2, signdpet == 3 ~ 1, signdpet == 4 ~ 0, TRUE ~ NA_real_),
      recod_avoidbuy = case_when(avoidbuy == 1 ~ 3, avoidbuy == 2 ~ 2, avoidbuy == 3 ~ 1, avoidbuy == 4 ~ 0, TRUE ~ NA_real_),
      recod_joindem  = case_when(joindem == 1 ~ 3, joindem == 2 ~ 2, joindem == 3 ~ 1, joindem == 4 ~ 0, TRUE ~ NA_real_),
      recod_attrally = case_when(attrally == 1 ~ 3, attrally == 2 ~ 2, attrally == 3 ~ 1, attrally == 4 ~ 0, TRUE ~ NA_real_),
      recod_cntctgov = case_when(cntctgov == 1 ~ 3, cntctgov == 2 ~ 2, cntctgov == 3 ~ 1, cntctgov == 4 ~ 0, TRUE ~ NA_real_),
      recod_polfunds = case_when(polfunds == 1 ~ 3, polfunds == 2 ~ 2, polfunds == 3 ~ 1, polfunds == 4 ~ 0, TRUE ~ NA_real_),
      recod_usemedia = case_when(usemedia == 1 ~ 3, usemedia == 2 ~ 2, usemedia == 3 ~ 1, usemedia == 4 ~ 0, TRUE ~ NA_real_),
      recod_interpol = case_when(interpol == 1 ~ 3, interpol == 2 ~ 2, interpol == 3 ~ 1, interpol == 4 ~ 0, TRUE ~ NA_real_),
      # Para actlaw (1=muy probable -> 3, 4=nada probable -> 0)
      recod_actlaw   = case_when(actlaw == 1 ~ 3, actlaw == 2 ~ 2, actlaw == 3 ~ 1, actlaw == 4 ~ 0, TRUE ~ NA_real_)
    ) %>%
    # Segundo, sumar las variables recodificadas y normalizar
    mutate(
      # Suma de los scores (rango 0-27). Si hay algún NA, el resultado será NA.
      score_bruto_normativo = rowSums(select(., starts_with("recod_"))),
      
      # Normalizar el score para que esté entre 0 y 1
      propensity_score = score_bruto_normativo / 27
    )
  
  # --- 4. Asignar el Score Calculado como un Nuevo Atributo a la Red ---
  set.vertex.attribute(gss_network, "propensity_score", df_propensity_attr$propensity_score)
  
  # --- 5. Guardar la Red Modificada ---
  output_filename <- file.path(networks_dir, sprintf("%s%03d.rds", output_prefix, i))
  saveRDS(gss_network, output_filename)
}

cat("Proceso completado. Todas las redes han sido actualizadas con el 'propensity_score'.\n\n")


# --- 6. Análisis y Visualización de UNA de las Redes Modificadas ---
cat("Generando gráficos de distribución para la primera red...\n")

# Cargar la primera red modificada para el análisis
first_mur_net_file <- file.path(networks_dir, sprintf("%s001.rds", output_prefix))
if (file.exists(first_mur_net_file)) {
  gss_network_mur <- readRDS(first_mur_net_file)
  
  # Extraer los atributos a un dataframe
  df_attributes <- data.frame(.ID = 1:network.size(gss_network_mur))
  for (attr in list.vertex.attributes(gss_network_mur)) {
    df_attributes[[attr]] <- get.vertex.attribute(gss_network_mur, attr)
  }
  
  # --- 6.1 Gráficos de Distribución de las Variables Originales ---
  plots_propensity_original <- list()
  for (p_var in propensity_ingredient_vars) {
    # Asegurar que la variable es un factor para el gráfico de barras
    df_attributes[[p_var]] <- factor(df_attributes[[p_var]], levels = 1:4)
    
    p <- ggplot(df_attributes, aes(x = .data[[p_var]], fill = .data[[p_var]])) +
      geom_bar() +
      scale_x_discrete(drop = FALSE) + # No eliminar categorías sin observaciones
      labs(title = paste("Var:", p_var), x = "Respuesta", y = "Frecuencia") +
      theme_minimal() +
      theme(legend.position = "none")
    plots_propensity_original[[p_var]] <- p
  }
  
  # Guardar los gráficos en un PDF
  pdf("trabajo_1_plots/gss_propensity_vars_distribution.pdf", width = 12, height = 9)
  do.call(grid.arrange, c(plots_propensity_original, ncol = 3))
  dev.off()
  
  # --- 6.2 Gráfico de Distribución del Score Final ---
  hist_score <- ggplot(df_attributes, aes(x = propensity_score)) +
    geom_histogram(bins = 28, fill = "skyblue", color = "black") + # 28 bins para 27 posibles valores + 0
    labs(title = "Distribución del Score de Propensión a la Acción Colectiva",
         x = "Propensity Score (0 a 1)", 
         y = "Frecuencia") +
    theme_minimal()
  
  print(hist_score)
  ggsave("trabajo_1_plots/gss_propensity_score_distribution.pdf", plot = hist_score, width = 8, height = 6)
  
  cat("Gráficos de distribución guardados en la carpeta 'trabajo_1_plots'.\n")
  
} else {
  warning("No se encontró la primera red modificada para generar los gráficos.")
}