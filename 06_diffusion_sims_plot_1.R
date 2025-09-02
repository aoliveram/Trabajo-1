# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script para Analizar Resultados y Generar CUATRO TIPOS de Heatmaps:
# 1. "Proporción de TODAS las Runs que son Exitosas Y Tienen Transición"
# 2. "Proporción Final Media de Adoptadores"
# 3. "Proporción Media de Nuevos Adoptantes por Elección Racional"
# 4. "Proporción Media de Nuevos Adoptantes por Influencia Social"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(viridis)

# --- Parámetros de Análisis ---
RESULTS_DIR <- "trabajo_1_files/diffusion_simulation_files/"
PLOTS_DIR <- "trabajo_1_plots/diffusion_simulation_plots/"
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP <- seq(0/12, 12/12, by = 1/12)
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6)
TAU_NORMAL_SD_SWEEP_LIST <- c(0.08, 0.12, 0.16, 0.20)
# NUM_SEED_RUNS_TOTAL se determinará a partir de los datos cargados

PHASE_TRANSITION_THRESHOLD_JUMP <- 0.50
SUCCESSFUL_DIFFUSION_THRESHOLD_PROP <- 0.50
MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL <- 0.80 # Para el heatmap de transiciones

cat("Cargando resultados crudos guardados (versión _2 con tipos de adopción)...\n")
# Asegúrate de cargar el archivo correcto que contiene las nuevas columnas
grand_raw_results_path <- paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_raw_results_all_sds_means.rds") # ASUMIENDO _2
if (!file.exists(grand_raw_results_path)) {
  stop("Archivo de resultados crudos no encontrado: ", grand_raw_results_path)
}
all_sds_raw_results_list_v2 <- readRDS(grand_raw_results_path)
cat("Resultados cargados.\n")

# --- Almacenamiento de Datos para Heatmaps ---
all_sds_transition_metric_heatmap_data_list <- list()
all_sds_avg_adoption_heatmap_data_list <- list()
all_sds_prop_rational_heatmap_data_list <- list() # NUEVO
all_sds_prop_social_heatmap_data_list <- list()   # NUEVO

# --- Procesamiento de Datos ---
for (current_tau_sd in TAU_NORMAL_SD_SWEEP_LIST) {
  sd_label <- paste0("sd_", sprintf("%.2f", current_tau_sd))
  cat(paste0("\nProcesando para Tau SD = ", current_tau_sd, "...\n"))
  
  current_sd_all_means_raw_results <- all_sds_raw_results_list_v2[[sd_label]]
  if (is.null(current_sd_all_means_raw_results)) {
    cat(paste0("  No hay datos para SD = ", current_tau_sd, ". Saltando.\n"))
    next
  }
  
  heatmap_data_transition_metric_this_sd_list <- list()
  heatmap_data_avg_adoption_this_sd_list <- list()
  heatmap_data_prop_rational_this_sd_list <- list() # NUEVO
  heatmap_data_prop_social_this_sd_list <- list()   # NUEVO
  
  for (current_threshold_mean in THRESHOLD_MEAN_SWEEP_LIST) {
    mean_label <- paste0("mean_", sprintf("%.2f", current_threshold_mean))
    cat(paste0("  Procesando para Mean τ = ", current_threshold_mean, " (dentro de SD = ", current_tau_sd, ")...\n"))
    
    raw_df_this_mean_sd <- current_sd_all_means_raw_results[[mean_label]]
    if (is.null(raw_df_this_mean_sd) || nrow(raw_df_this_mean_sd) == 0) {
      cat(paste0("    No hay datos para Mean τ = ", current_threshold_mean, ". Saltando.\n"))
      next
    }
    
    NUM_RUNS_THIS_COMBO_ACTUAL <- length(unique(raw_df_this_mean_sd$run_id))
    if(NUM_RUNS_THIS_COMBO_ACTUAL == 0) {
      cat(paste0("    DataFrame raw_df_this_mean_sd vacío para Mean τ = ", current_threshold_mean, ". Saltando.\n"))
      next
    }
    
    # --- 1. Pre-cálculo para el heatmap de TRANSICIONES ---
    run_summary_info_for_transitions <- raw_df_this_mean_sd %>%
      group_by(run_id, social_distance_h, innovation_iul_Gamma) %>%
      summarise(
        adopters_prop_at_cell = num_adopters / N_nodes_actual,
        # Guardar las columnas necesarias para los nuevos heatmaps
        num_adopted_rational_at_cell = first(num_adopted_rational), # Debería ser constante para esta celda (Γ,h) y run_id
        num_adopted_social_at_cell = first(num_adopted_social),
        num_adopters_at_cell = first(num_adopters),
        initial_cluster_size_at_cell = first(initial_cluster_size),
        n_nodes_at_cell = first(N_nodes_actual),
        .groups = 'drop'
      ) %>%
      mutate(
        is_successful_at_cell = adopters_prop_at_cell >= SUCCESSFUL_DIFFUSION_THRESHOLD_PROP
      ) %>%
      arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
      group_by(run_id, social_distance_h) %>%
      mutate(
        jump_at_step = adopters_prop_at_cell - lag(adopters_prop_at_cell, default = 0),
        is_transition_at_step = ifelse(jump_at_step >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
      ) %>%
      group_by(run_id, social_distance_h) %>%
      mutate(
        first_transition_IUL_for_series = if (any(is_transition_at_step == 1, na.rm=TRUE)) {
          min(innovation_iul_Gamma[which(is_transition_at_step == 1)])
        } else { NA_real_ }
      ) %>%
      ungroup()
    
    # --- 2. Pre-cálculo para el heatmap de ADOPCIÓN MEDIA (y nuevas proporciones) ---
    # Esta agregación ahora también calculará las proporciones medias de adopción racional/social
    aggregated_data_per_cell_df <- raw_df_this_mean_sd %>%
      mutate(
        new_adopters = num_adopters - initial_cluster_size,
        prop_rational_of_new = ifelse(new_adopters > 0, num_adopted_rational / new_adopters, 0), # Evitar NaN si new_adopters es 0
        prop_social_of_new = ifelse(new_adopters > 0, num_adopted_social / new_adopters, 0)
      ) %>%
      group_by(innovation_iul_Gamma, social_distance_h) %>%
      summarise(
        mean_final_adopters_prop = mean(num_adopters / N_nodes_actual, na.rm = TRUE),
        mean_prop_rational_of_new = mean(prop_rational_of_new, na.rm = TRUE), # NUEVO
        mean_prop_social_of_new = mean(prop_social_of_new, na.rm = TRUE),     # NUEVO
        .groups = 'drop'
      )
    # Corregir NaNs si todos los new_adopters fueron 0 para una celda, resultando en mean(numeric(0))
    aggregated_data_per_cell_df$mean_prop_rational_of_new[is.nan(aggregated_data_per_cell_df$mean_prop_rational_of_new)] <- 0
    aggregated_data_per_cell_df$mean_prop_social_of_new[is.nan(aggregated_data_per_cell_df$mean_prop_social_of_new)] <- 0
    
    
    # --- Llenar datos para los 4 heatmaps ---
    panel_heatmap_data_transition_metric_list <- list()
    panel_heatmap_data_avg_adoption_list <- list()
    panel_heatmap_data_prop_rational_list <- list() # NUEVO
    panel_heatmap_data_prop_social_list <- list()   # NUEVO
    
    for (iul_val in IUL_VALUES_SWEEP) {
      for (h_val in H_VALUES_SWEEP) {
        
        # --- Para Heatmap de Transiciones (como antes) ---
        runs_data_for_this_cell_transitions <- run_summary_info_for_transitions %>%
          filter(innovation_iul_Gamma == iul_val, social_distance_h == h_val)
        prop_transition_metric_this_cell <- NA_real_
        if (nrow(runs_data_for_this_cell_transitions) > 0) {
          num_successful_runs_in_this_cell_t <- sum(runs_data_for_this_cell_transitions$is_successful_at_cell, na.rm = TRUE)
          prop_successful_runs_in_this_cell_t <- num_successful_runs_in_this_cell_t / NUM_RUNS_THIS_COMBO_ACTUAL
          if (prop_successful_runs_in_this_cell_t >= MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL) {
            runs_successful_and_transitioned_for_cell <- runs_data_for_this_cell_transitions %>%
              filter(is_successful_at_cell == TRUE & !is.na(first_transition_IUL_for_series) &
                       first_transition_IUL_for_series <= iul_val)
            prop_transition_metric_this_cell <- nrow(runs_successful_and_transitioned_for_cell) / NUM_RUNS_THIS_COMBO_ACTUAL
          }
        }
        panel_heatmap_data_transition_metric_list[[length(panel_heatmap_data_transition_metric_list) + 1]] <- data.frame(
          innovation_iul_Gamma = iul_val, social_distance_h = h_val,
          proportion_value_to_plot = prop_transition_metric_this_cell, 
          tau_mean_param = current_threshold_mean
        )
        
        # --- Para Heatmap de Adopción Media (como antes) ---
        current_cell_avg_adoption_df <- aggregated_data_per_cell_df %>%
          filter(innovation_iul_Gamma == iul_val, social_distance_h == h_val)
        mean_adoption_to_plot_this_cell <- NA_real_
        if(nrow(current_cell_avg_adoption_df) > 0){
          mean_adoption_to_plot_this_cell <- current_cell_avg_adoption_df$mean_final_adopters_prop
        }
        panel_heatmap_data_avg_adoption_list[[length(panel_heatmap_data_avg_adoption_list) + 1]] <- data.frame(
          innovation_iul_Gamma = iul_val, social_distance_h = h_val,
          mean_adopters_prop_to_plot = mean_adoption_to_plot_this_cell,
          tau_mean_param = current_threshold_mean
        )
        
        # --- NUEVO: Para Heatmap de Proporción Racional ---
        mean_prop_rational_to_plot_this_cell <- NA_real_
        if(nrow(current_cell_avg_adoption_df) > 0){ # Reutilizamos current_cell_avg_adoption_df
          mean_prop_rational_to_plot_this_cell <- current_cell_avg_adoption_df$mean_prop_rational_of_new
        }
        panel_heatmap_data_prop_rational_list[[length(panel_heatmap_data_prop_rational_list) + 1]] <- data.frame(
          innovation_iul_Gamma = iul_val, social_distance_h = h_val,
          mean_prop_rational_to_plot = mean_prop_rational_to_plot_this_cell,
          tau_mean_param = current_threshold_mean
        )
        
        # --- NUEVO: Para Heatmap de Proporción Social ---
        mean_prop_social_to_plot_this_cell <- NA_real_
        if(nrow(current_cell_avg_adoption_df) > 0){ # Reutilizamos current_cell_avg_adoption_df
          mean_prop_social_to_plot_this_cell <- current_cell_avg_adoption_df$mean_prop_social_of_new
        }
        panel_heatmap_data_prop_social_list[[length(panel_heatmap_data_prop_social_list) + 1]] <- data.frame(
          innovation_iul_Gamma = iul_val, social_distance_h = h_val,
          mean_prop_social_to_plot = mean_prop_social_to_plot_this_cell,
          tau_mean_param = current_threshold_mean
        )
      }
    }
    heatmap_data_transition_metric_this_sd_list[[mean_label]] <- bind_rows(panel_heatmap_data_transition_metric_list)
    heatmap_data_avg_adoption_this_sd_list[[mean_label]] <- bind_rows(panel_heatmap_data_avg_adoption_list)
    heatmap_data_prop_rational_this_sd_list[[mean_label]] <- bind_rows(panel_heatmap_data_prop_rational_list) # NUEVO
    heatmap_data_prop_social_this_sd_list[[mean_label]] <- bind_rows(panel_heatmap_data_prop_social_list)   # NUEVO
    cat(paste0("    Datos para los cuatro heatmaps procesados para Mean τ = ", current_threshold_mean, ".\n"))
  } 
  
  # --- Guardar datos pre-plot para el CURRENT SD ---
  if (length(heatmap_data_transition_metric_this_sd_list) > 0) { # Transiciones
    valid_keys_t <- !sapply(heatmap_data_transition_metric_this_sd_list, is.null)
    if(sum(valid_keys_t) > 0) {
      df_heatmap_this_sd_transitions <- bind_rows(heatmap_data_transition_metric_this_sd_list[valid_keys_t])
      df_heatmap_this_sd_transitions$tau_sd_param <- current_tau_sd 
      all_sds_transition_metric_heatmap_data_list[[sd_label]] <- df_heatmap_this_sd_transitions 
    }
  }
  if (length(heatmap_data_avg_adoption_this_sd_list) > 0) { # Adopción Media
    valid_keys_a <- !sapply(heatmap_data_avg_adoption_this_sd_list, is.null)
    if(sum(valid_keys_a) > 0) {
      df_heatmap_this_sd_avg_adoption <- bind_rows(heatmap_data_avg_adoption_this_sd_list[valid_keys_a])
      df_heatmap_this_sd_avg_adoption$tau_sd_param <- current_tau_sd
      all_sds_avg_adoption_heatmap_data_list[[sd_label]] <- df_heatmap_this_sd_avg_adoption
    }
  }
  if (length(heatmap_data_prop_rational_this_sd_list) > 0) { # Prop. Racional
    valid_keys_r <- !sapply(heatmap_data_prop_rational_this_sd_list, is.null)
    if(sum(valid_keys_r) > 0) {
      df_heatmap_this_sd_prop_rational <- bind_rows(heatmap_data_prop_rational_this_sd_list[valid_keys_r])
      df_heatmap_this_sd_prop_rational$tau_sd_param <- current_tau_sd
      all_sds_prop_rational_heatmap_data_list[[sd_label]] <- df_heatmap_this_sd_prop_rational
    }
  }
  if (length(heatmap_data_prop_social_this_sd_list) > 0) { # Prop. Social
    valid_keys_s <- !sapply(heatmap_data_prop_social_this_sd_list, is.null)
    if(sum(valid_keys_s) > 0) {
      df_heatmap_this_sd_prop_social <- bind_rows(heatmap_data_prop_social_this_sd_list[valid_keys_s])
      df_heatmap_this_sd_prop_social$tau_sd_param <- current_tau_sd
      all_sds_prop_social_heatmap_data_list[[sd_label]] <- df_heatmap_this_sd_prop_social
    }
  }
} 

cat("\nProcesamiento de datos para todos los SDs completado.\n")

# --- Generación de PLOTS (en bucles separados por tipo de plot) ---

# Función auxiliar para crear un plot de heatmap (para reducir duplicación)
create_heatmap_plot <- function(df_plot_data, fill_variable_name, legend_title, plot_main_title_prefix, current_sd_val, num_runs_val, facet_var_name, color_option="plasma", subtitle_note="") {
  if (nrow(df_plot_data) == 0 || !any(!is.na(df_plot_data[[fill_variable_name]]))) {
    cat(paste0("  Skipping plot for SD = ", current_sd_val, " (", plot_main_title_prefix, ") due to no plottable data.\n"))
    return(NULL)
  }
  
  df_plot_data$social_distance_h_factor <- factor(sprintf("%.2f", df_plot_data$social_distance_h))
  df_plot_data$facet_label <- factor(paste0("Mean τ = ", sprintf("%.2f", df_plot_data[[facet_var_name]])))
  
  plot_obj <- ggplot(df_plot_data, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = .data[[fill_variable_name]])) +
    geom_tile(color = "grey80", lwd = 0.1) + 
    scale_fill_viridis_c(name = legend_title, limits = c(0, 1), option=color_option, n.breaks=6, na.value = "white") +
    facet_wrap(~facet_label, ncol = 2) + 
    labs(x = expression(paste("Intrinsic Utility Level - IUL (", Gamma, ")")),
         y = "Maximum Social Closeness - MSP (h)", 
         title = paste(plot_main_title_prefix, "- Tau SD =", sprintf("%.2f", current_sd_val)),
         subtitle = paste("Thresholds ~ N(μ=var, σ=", sprintf("%.2f", current_sd_val), "), ", num_runs_val, 
                          " runs. ", subtitle_note, sep="")
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0,0)) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=11), 
          plot.subtitle = element_text(hjust = 0.5, size=7),
          axis.text.x = element_text(angle = 45, hjust = 1, size=7), axis.text.y = element_text(size=7), 
          legend.position = "right", legend.title = element_text(size = 7), legend.text = element_text(size = 6),
          legend.key.size = unit(0.7, "lines"), strip.text = element_text(face="bold", size=8))
  return(plot_obj)
}


# 1. Plot para PHASE TRANSITION
cat("Generando plots para 'Transition Metric'...\n")
grand_combined_transition_heatmap_df <- bind_rows(all_sds_transition_metric_heatmap_data_list)
if(nrow(grand_combined_transition_heatmap_df) > 0){
  for (current_tau_sd_plot in TAU_NORMAL_SD_SWEEP_LIST) {
    df_plot_tm <- grand_combined_transition_heatmap_df %>% filter(tau_sd_param == current_tau_sd_plot)
    subtitle_tm <- paste0("Tiles shown if >= ", MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL*100,"% runs in cell successful.")
    plot_obj_tm <- create_heatmap_plot(df_plot_tm, "proportion_value_to_plot", 
                                       "Prop. of All Runs:\nSuccessful in Cell\n& Phase Transition",
                                       "Phase Transition Maps (Cell-Successful)", 
                                       current_tau_sd_plot, NUM_RUNS_THIS_COMBO_ACTUAL, "tau_mean_param", 
                                       "plasma", subtitle_tm)
    if(!is.null(plot_obj_tm)){
      filename_tm <- paste0(PLOTS_DIR, "phase_transition_sd", sprintf("%.2f", current_tau_sd_plot), "_means03-06.pdf")
      ggsave(filename_tm, plot_obj_tm, width = 6.5, height = 4.5)
      cat(paste0("  Saved plot (Transition Metric): ", filename_tm, "\n"))
    }
  }
} else {cat("No data for 'Transition Metric' heatmaps.\n")}

# 2. Plot para AVERAGE ADOPTION
cat("Generando plots para 'Average Adoption'...\n")
grand_combined_avg_adoption_heatmap_df <- bind_rows(all_sds_avg_adoption_heatmap_data_list)
if(nrow(grand_combined_avg_adoption_heatmap_df) > 0){
  for (current_tau_sd_plot in TAU_NORMAL_SD_SWEEP_LIST) {
    df_plot_aa <- grand_combined_avg_adoption_heatmap_df %>% filter(tau_sd_param == current_tau_sd_plot)
    plot_obj_aa <- create_heatmap_plot(df_plot_aa, "mean_adopters_prop_to_plot",
                                       "Mean Final\nAdopter Prop.",
                                       "Mean Final Adopter Proportion Maps",
                                       current_tau_sd_plot, NUM_RUNS_THIS_COMBO_ACTUAL, "tau_mean_param",
                                       "cividis") # No subtitle note needed for this one
    if(!is.null(plot_obj_aa)){
      filename_aa <- paste0(PLOTS_DIR, "avg_adoption_sd", sprintf("%.2f", current_tau_sd_plot), "_means03-06.pdf")
      ggsave(filename_aa, plot_obj_aa, width = 6.5, height = 4.5)
      cat(paste0("  Saved plot (Average Adoption): ", filename_aa, "\n"))
    }
  }
} else {cat("No data for 'Average Adoption' heatmaps.\n")}

# 3. Plot para PROP. RACIONAL
cat("Generando plots para 'Proportion Rational Adoption'...\n")
grand_combined_prop_rational_df <- bind_rows(all_sds_prop_rational_heatmap_data_list)
if(nrow(grand_combined_prop_rational_df) > 0){
  for (current_tau_sd_plot in TAU_NORMAL_SD_SWEEP_LIST) {
    df_plot_pr <- grand_combined_prop_rational_df %>% filter(tau_sd_param == current_tau_sd_plot)
    plot_obj_pr <- create_heatmap_plot(df_plot_pr, "mean_prop_rational_to_plot",
                                       "Mean Prop. New Adopters\nvia Rational Choice",
                                       "Proportion Rational Adoption Maps",
                                       current_tau_sd_plot, NUM_RUNS_THIS_COMBO_ACTUAL, "tau_mean_param",
                                       "magma") # Another palette
    if(!is.null(plot_obj_pr)){
      filename_pr <- paste0(PLOTS_DIR, "avg_adoption_rational_sd", sprintf("%.2f", current_tau_sd_plot), "_means03-06.pdf")
      ggsave(filename_pr, plot_obj_pr, width = 6.5, height = 4.5)
      cat(paste0("  Saved plot (Prop. Rational): ", filename_pr, "\n"))
    }
  }
} else {cat("No data for 'Proportion Rational Adoption' heatmaps.\n")}

# 4. Plot para PROP. SOCIAL
cat("Generando plots para 'Proportion Social Adoption'...\n")
grand_combined_prop_social_df <- bind_rows(all_sds_prop_social_heatmap_data_list)
if(nrow(grand_combined_prop_social_df) > 0){
  for (current_tau_sd_plot in TAU_NORMAL_SD_SWEEP_LIST) {
    df_plot_ps <- grand_combined_prop_social_df %>% filter(tau_sd_param == current_tau_sd_plot)
    plot_obj_ps <- create_heatmap_plot(df_plot_ps, "mean_prop_social_to_plot",
                                       "Mean Prop. New Adopters\nvia Social Influence",
                                       "Proportion Social Influence Maps",
                                       current_tau_sd_plot, NUM_RUNS_THIS_COMBO_ACTUAL, "tau_mean_param",
                                       "inferno") # Yet another palette
    if(!is.null(plot_obj_ps)){
      filename_ps <- paste0(PLOTS_DIR, "avg_adoption_socInfluenc_sd", sprintf("%.2f", current_tau_sd_plot), "_means03-06.pdf")
      ggsave(filename_ps, plot_obj_ps, width = 6.5, height = 4.5)
      cat(paste0("  Saved plot (Prop. Social): ", filename_ps, "\n"))
    }
  }
} else {cat("No data for 'Proportion Social Adoption' heatmaps.\n")}


# --- Guardar los datos pre-plot finales ---
saveRDS(all_sds_transition_metric_heatmap_data_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_FINAL_METRIC_heatmap_data_all_sds_v2.rds"))
saveRDS(all_sds_avg_adoption_heatmap_data_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_AVG_ADOPTION_heatmap_data_all_sds_v2.rds"))
saveRDS(all_sds_prop_rational_heatmap_data_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_PROP_RATIONAL_heatmap_data_all_sds_v2.rds"))
saveRDS(all_sds_prop_social_heatmap_data_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_PROP_SOCIAL_heatmap_data_all_sds_v2.rds"))

cat("Todos los análisis y guardados completados.\n")