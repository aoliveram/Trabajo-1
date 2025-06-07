# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script para Analizar Resultados y Generar DOS TIPOS de Heatmaps:
# 1. "Proporción de TODAS las Runs que son Exitosas Y Tienen Transición"
# 2. "Proporción Final Media de Adoptadores"
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
# NUM_SEED_RUNS_TOTAL <- 100 # Lo tomaré del número de run_id únicos por (mean,sd)
# Si en tus datos hay 96 runs, eso se usará.

PHASE_TRANSITION_THRESHOLD_JUMP <- 0.50
SUCCESSFUL_DIFFUSION_THRESHOLD_PROP <- 0.50
MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL <- 0.80 # Para el heatmap de transiciones

cat("Cargando resultados crudos guardados...\n")
grand_raw_results_path <- paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_raw_results_all_sds_means.rds")
if (!file.exists(grand_raw_results_path)) {
  stop("Archivo de resultados crudos no encontrado: ", grand_raw_results_path)
}
all_sds_raw_results_list <- readRDS(grand_raw_results_path)
cat("Resultados cargados.\n")

# --- Almacenamiento de Datos para Heatmaps ---
all_sds_transition_metric_heatmap_data_list <- list() # Para el heatmap de transiciones
all_sds_avg_adoption_heatmap_data_list <- list()    # Para el heatmap de adopción media

# --- Procesamiento de Datos ---
for (current_tau_sd in TAU_NORMAL_SD_SWEEP_LIST) {
  sd_label <- paste0("sd_", sprintf("%.2f", current_tau_sd))
  cat(paste0("\nProcesando para Tau SD = ", current_tau_sd, "...\n"))
  
  current_sd_all_means_raw_results <- all_sds_raw_results_list[[sd_label]]
  if (is.null(current_sd_all_means_raw_results)) {
    cat(paste0("  No hay datos para SD = ", current_tau_sd, ". Saltando.\n"))
    next
  }
  
  heatmap_data_transition_metric_this_sd_list <- list()
  heatmap_data_avg_adoption_this_sd_list <- list()
  
  for (current_threshold_mean in THRESHOLD_MEAN_SWEEP_LIST) {
    mean_label <- paste0("mean_", sprintf("%.2f", current_threshold_mean))
    cat(paste0("  Procesando para Mean τ = ", current_threshold_mean, " (dentro de SD = ", current_tau_sd, ")...\n"))
    
    raw_df_this_mean_sd <- current_sd_all_means_raw_results[[mean_label]]
    if (is.null(raw_df_this_mean_sd) || nrow(raw_df_this_mean_sd) == 0) {
      cat(paste0("    No hay datos para Mean τ = ", current_threshold_mean, ". Saltando.\n"))
      next
    }
    
    # Determinar el número real de runs para esta combinación (en caso de que no sean 100)
    ACTUAL_NUM_RUNS_THIS_COMBO <- length(unique(raw_df_this_mean_sd$run_id))
    if(ACTUAL_NUM_RUNS_THIS_COMBO == 0) { # Si el dataframe está vacío después de todo
      cat(paste0("    DataFrame raw_df_this_mean_sd vacío para Mean τ = ", current_threshold_mean, ". Saltando.\n"))
      next
    }
    
    
    # --- 1. Pre-cálculo para el heatmap de TRANSICIONES ---
    run_summary_info_for_transitions <- raw_df_this_mean_sd %>%
      group_by(run_id, social_distance_h, innovation_iul_Gamma) %>%
      summarise(
        adopters_prop_at_cell = num_adopters / N_nodes_actual,
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
      group_by(run_id, social_distance_h) %>% # Re-agrupar para el siguiente mutate
      mutate(
        first_transition_IUL_for_series = if (any(is_transition_at_step == 1, na.rm=TRUE)) {
          min(innovation_iul_Gamma[which(is_transition_at_step == 1)])
        } else { NA_real_ }
      ) %>%
      ungroup()
    
    # --- 2. Pre-cálculo para el heatmap de ADOPCIÓN MEDIA ---
    avg_adoption_per_cell_df <- raw_df_this_mean_sd %>%
      group_by(innovation_iul_Gamma, social_distance_h) %>%
      summarise(
        mean_final_adopters_prop = mean(num_adopters / N_nodes_actual, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # --- Llenar datos para ambos heatmaps ---
    panel_heatmap_data_transition_metric_list <- list()
    panel_heatmap_data_avg_adoption_list <- list()
    
    for (iul_val in IUL_VALUES_SWEEP) {
      for (h_val in H_VALUES_SWEEP) {
        
        # --- Para Heatmap de Transiciones ---
        runs_data_for_this_cell_transitions <- run_summary_info_for_transitions %>%
          filter(innovation_iul_Gamma == iul_val, social_distance_h == h_val)
        
        prop_transition_metric_this_cell <- NA_real_
        
        if (nrow(runs_data_for_this_cell_transitions) > 0) {
          num_successful_runs_in_this_cell <- sum(runs_data_for_this_cell_transitions$is_successful_at_cell, na.rm = TRUE)
          prop_successful_runs_in_this_cell <- num_successful_runs_in_this_cell / ACTUAL_NUM_RUNS_THIS_COMBO
          
          if (prop_successful_runs_in_this_cell >= MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL) {
            runs_successful_and_transitioned_for_cell <- runs_data_for_this_cell_transitions %>%
              filter(is_successful_at_cell == TRUE &
                       !is.na(first_transition_IUL_for_series) &
                       first_transition_IUL_for_series <= iul_val)
            
            num_successful_and_transitioned_in_cell <- nrow(runs_successful_and_transitioned_for_cell)
            prop_transition_metric_this_cell <- num_successful_and_transitioned_in_cell / ACTUAL_NUM_RUNS_THIS_COMBO
          }
        }
        panel_heatmap_data_transition_metric_list[[length(panel_heatmap_data_transition_metric_list) + 1]] <- data.frame(
          innovation_iul_Gamma = iul_val, social_distance_h = h_val,
          proportion_value_to_plot = prop_transition_metric_this_cell, 
          tau_mean_param = current_threshold_mean
        )
        
        # --- Para Heatmap de Adopción Media ---
        avg_adoption_this_cell_df <- avg_adoption_per_cell_df %>%
          filter(innovation_iul_Gamma == iul_val, social_distance_h == h_val)
        
        mean_adoption_to_plot_this_cell <- NA_real_
        if(nrow(avg_adoption_this_cell_df) > 0){
          mean_adoption_to_plot_this_cell <- avg_adoption_this_cell_df$mean_final_adopters_prop
        }
        panel_heatmap_data_avg_adoption_list[[length(panel_heatmap_data_avg_adoption_list) + 1]] <- data.frame(
          innovation_iul_Gamma = iul_val, social_distance_h = h_val,
          mean_adopters_prop_to_plot = mean_adoption_to_plot_this_cell,
          tau_mean_param = current_threshold_mean
        )
      }
    }
    heatmap_data_transition_metric_this_sd_list[[mean_label]] <- bind_rows(panel_heatmap_data_transition_metric_list)
    heatmap_data_avg_adoption_this_sd_list[[mean_label]] <- bind_rows(panel_heatmap_data_avg_adoption_list)
    cat(paste0("    Datos para ambos heatmaps procesados para Mean τ = ", current_threshold_mean, ".\n"))
  } 
  
  # --- Guardar datos para el heatmap de TRANSICIONES para el CURRENT SD ---
  if (length(heatmap_data_transition_metric_this_sd_list) > 0) {
    valid_keys_t <- !sapply(heatmap_data_transition_metric_this_sd_list, is.null)
    if(sum(valid_keys_t) > 0) {
      df_heatmap_this_sd_transitions <- bind_rows(heatmap_data_transition_metric_this_sd_list[valid_keys_t])
      df_heatmap_this_sd_transitions$tau_sd_param <- current_tau_sd 
      all_sds_transition_metric_heatmap_data_list[[sd_label]] <- df_heatmap_this_sd_transitions 
    }
  }
  # --- Guardar datos para el heatmap de ADOPCIÓN MEDIA para el CURRENT SD ---
  if (length(heatmap_data_avg_adoption_this_sd_list) > 0) {
    valid_keys_a <- !sapply(heatmap_data_avg_adoption_this_sd_list, is.null)
    if(sum(valid_keys_a) > 0) {
      df_heatmap_this_sd_avg_adoption <- bind_rows(heatmap_data_avg_adoption_this_sd_list[valid_keys_a])
      df_heatmap_this_sd_avg_adoption$tau_sd_param <- current_tau_sd
      all_sds_avg_adoption_heatmap_data_list[[sd_label]] <- df_heatmap_this_sd_avg_adoption
    }
  }
} # Fin bucle TAU_NORMAL_SD_SWEEP_LIST

cat("\nProcesamiento de datos para todos los SDs completado.\n")

# --- Generación de PLOTS ---

# 1. Plot para TRANSITION METRIC
cat("Generando plots para 'Transition Metric'...\n")
grand_combined_transition_heatmap_df <- bind_rows(all_sds_transition_metric_heatmap_data_list)

if (nrow(grand_combined_transition_heatmap_df) > 0) {
  for (current_tau_sd_plot in TAU_NORMAL_SD_SWEEP_LIST) {
    df_plot_tm <- grand_combined_transition_heatmap_df %>% filter(tau_sd_param == current_tau_sd_plot)
    
    if (nrow(df_plot_tm) > 0 && any(!is.na(df_plot_tm$proportion_value_to_plot))) {
      df_plot_tm$social_distance_h_factor <- factor(sprintf("%.2f", df_plot_tm$social_distance_h))
      df_plot_tm$tau_mean_facet_label <- factor(paste0("Mean τ = ", sprintf("%.2f", df_plot_tm$tau_mean_param)))
      
      legend_title_tm <- "Proportion of\nruns with\nPhase Transition"
      plot_title_tm <- paste("Phase Transition Maps (Cell-Successful) - Tau SD =", sprintf("%.2f", current_tau_sd_plot))
      plot_subtitle_tm <- paste("Thresholds ~ N(mu=var, sigma=", sprintf("%.2f", current_tau_sd_plot), "), ", ACTUAL_NUM_RUNS_THIS_COMBO, # Usar el último valor de ACTUAL_NUM_RUNS_THIS_COMBO o el NUM_SEED_RUNS_TOTAL global
                                " runs. Tiles shown if >= ", MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL*100,"% runs in cell successful.", sep="")
      
      
      plot_obj_tm <- ggplot(df_plot_tm, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = proportion_value_to_plot)) +
        geom_tile(color = "grey80", lwd = 0.1) + 
        scale_fill_viridis_c(name = legend_title_tm, limits = c(0, 1), option="plasma", n.breaks=6, na.value = "white") +
        facet_wrap(~tau_mean_facet_label, ncol = 2) + 
        labs(x = expression(paste("Intrinsic Utility Level - IUL (", Gamma, ")")),
             y = "Maximum Social Closeness - MSP (h)", title = plot_title_tm, subtitle = plot_subtitle_tm) +
        scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0,0)) +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=11), plot.subtitle = element_text(hjust = 0.5, size=7),
              axis.text.x = element_text(angle = 45, hjust = 1, size=7), axis.text.y = element_text(size=7), 
              legend.position = "right", legend.title = element_text(size = 7), legend.text = element_text(size = 6),
              legend.key.size = unit(0.7, "lines"), strip.text = element_text(face="bold", size=8))
      
      filename_tm <- paste0(PLOTS_DIR, "phase_transition_sd", sprintf("%.2f", current_tau_sd_plot), "_means03-06.pdf")
      ggsave(filename_tm, plot_obj_tm, width = 6.5, height = 4.5)
      cat(paste0("  Saved plot (Transition Metric): ", filename_tm, "\n"))
    } else {
      cat(paste0("  Skipping plot (Transition Metric) for SD = ", current_tau_sd_plot, " due to no plottable data.\n"))
    }
  }
} else {
  cat("No data available for 'Transition Metric' heatmaps.\n")
}

# 2. Plot para AVERAGE ADOPTION
cat("Generando plots para 'Average Adoption'...\n")
grand_combined_avg_adoption_heatmap_df <- bind_rows(all_sds_avg_adoption_heatmap_data_list)

if(nrow(grand_combined_avg_adoption_heatmap_df) > 0){
  for (current_tau_sd_plot in TAU_NORMAL_SD_SWEEP_LIST) {
    df_plot_aa <- grand_combined_avg_adoption_heatmap_df %>% filter(tau_sd_param == current_tau_sd_plot)
    
    if (nrow(df_plot_aa) > 0 && any(!is.na(df_plot_aa$mean_adopters_prop_to_plot))) {
      df_plot_aa$social_distance_h_factor <- factor(sprintf("%.2f", df_plot_aa$social_distance_h))
      df_plot_aa$tau_mean_facet_label <- factor(paste0("Mean τ = ", sprintf("%.2f", df_plot_aa$tau_mean_param)))
      
      legend_title_aa <- "Mean Final\nAdopter Prop."
      plot_title_aa <- paste("Mean Final Adopter Proportion Maps - Tau SD =", sprintf("%.2f", current_tau_sd_plot))
      plot_subtitle_aa <- paste("Thresholds ~ N(mu=var, sigma=", sprintf("%.2f", current_tau_sd_plot), "), ", ACTUAL_NUM_RUNS_THIS_COMBO, " runs per (Gamma,h) per panel", sep="")
      
      
      plot_obj_aa <- ggplot(df_plot_aa, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = mean_adopters_prop_to_plot)) +
        geom_tile(color = "grey80", lwd = 0.1) +
        scale_fill_viridis_c(name = legend_title_aa, limits = c(0, 1), option="cividis", n.breaks=6, na.value = "white") + # Usar otra paleta para diferenciar
        facet_wrap(~tau_mean_facet_label, ncol = 2) +
        labs(x = expression(paste("Intrinsic Utility Level - IUL (", Gamma, ")")),
             y = "Maximum Social Closeness - MSP (h)", title = plot_title_aa, subtitle = plot_subtitle_aa) +
        scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0,0)) +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=11), plot.subtitle = element_text(hjust = 0.5, size=7),
              axis.text.x = element_text(angle = 45, hjust = 1, size=7), axis.text.y = element_text(size=7),
              legend.position = "right", legend.title = element_text(size = 7), legend.text = element_text(size = 6),
              legend.key.size = unit(0.7, "lines"), strip.text = element_text(face="bold", size=8))
      
      filename_aa <- paste0(PLOTS_DIR, "avg_adoption_sd", sprintf("%.2f", current_tau_sd_plot), "_means03-06.pdf")
      ggsave(filename_aa, plot_obj_aa, width = 6.5, height = 4.5)
      cat(paste0("  Saved plot (Average Adoption): ", filename_aa, "\n"))
    } else {
      cat(paste0("  Skipping plot (Average Adoption) for SD = ", current_tau_sd_plot, " due to no plottable data.\n"))
    }
  }
} else {
  cat("No data available for 'Average Adoption' heatmaps.\n")
}


# --- Guardar los datos pre-plot finales ---
saveRDS(all_sds_transition_metric_heatmap_data_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_FINAL_METRIC_heatmap_data_all_sds.rds"))
saveRDS(all_sds_avg_adoption_heatmap_data_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_AVG_ADOPTION_heatmap_data_all_sds.rds"))

cat("Todos los análisis y guardados completados.\n")