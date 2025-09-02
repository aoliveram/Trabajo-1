# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script para Analizar Resultados CRUDOS y Generar PDF Consolidados con 16 Heatmaps
# Cada PDF corresponde a un SD_tau.
# Dentro del PDF: Filas por MEAN_tau, Columnas por TIPO de Heatmap.
# PARTE DE: phase_transition_GRAND_COMBINED_raw_results_all_sds_means_2.rds
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
SUCCESSFUL_DIFFUSION_THRESHOLD_PROP <- 0.50 # Para definir éxito EN CELDA
MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL <- 0.80 # Para el heatmap de transiciones

cat("Cargando resultados crudos guardados (versión _2 con tipos de adopción)...\n")
grand_raw_results_path <- paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_raw_results_all_sds_means.rds")
if (!file.exists(grand_raw_results_path)) {
  stop("Archivo de resultados crudos no encontrado: ", grand_raw_results_path)
}
all_sds_raw_results_list_v2 <- readRDS(grand_raw_results_path) # Esta es la lista de listas de dataframes raw
cat("Resultados crudos cargados.\n")


# --- PASO 0: Procesar los datos crudos para generar los 4 dataframes de heatmap ---
# Estas listas almacenarán los dataframes listos para plotear, un elemento por SD_tau
# Cada elemento será un dataframe que contiene datos para todas las MEDIAS de tau para ese SD_tau
all_sds_transition_metric_heatmap_df_list <- list()
all_sds_avg_adoption_heatmap_df_list    <- list()
all_sds_prop_rational_heatmap_df_list   <- list()
all_sds_prop_social_heatmap_df_list     <- list()

cat("\nPre-procesando datos crudos para todos los heatmaps...\n")

for (current_tau_sd_proc in TAU_NORMAL_SD_SWEEP_LIST) {
  sd_label_proc <- paste0("sd_", sprintf("%.2f", current_tau_sd_proc))
  cat(paste0("  Pre-procesando datos para Tau SD = ", current_tau_sd_proc, "...\n"))
  
  current_sd_all_means_raw_results_proc <- all_sds_raw_results_list_v2[[sd_label_proc]]
  if (is.null(current_sd_all_means_raw_results_proc)) {
    cat(paste0("    No hay datos crudos para SD = ", current_tau_sd_proc, ". Saltando pre-procesamiento.\n"))
    next
  }
  
  # Listas para almacenar los datos de heatmap de cada media DENTRO de este SD
  heatmap_data_tm_this_sd_list_proc <- list()
  heatmap_data_aa_this_sd_list_proc <- list()
  heatmap_data_pr_this_sd_list_proc <- list()
  heatmap_data_ps_this_sd_list_proc <- list()
  
  for (current_threshold_mean_proc in THRESHOLD_MEAN_SWEEP_LIST) {
    mean_label_proc <- paste0("mean_", sprintf("%.2f", current_threshold_mean_proc))
    # cat(paste0("    Pre-procesando para Mean τ = ", current_threshold_mean_proc, "...\n")) # Puede ser muy verboso
    
    raw_df_this_mean_sd_proc <- current_sd_all_means_raw_results_proc[[mean_label_proc]]
    if (is.null(raw_df_this_mean_sd_proc) || nrow(raw_df_this_mean_sd_proc) == 0) {
      next # Saltar si no hay datos para esta media
    }
    
    NUM_RUNS_THIS_COMBO_ACTUAL_PROC <- length(unique(raw_df_this_mean_sd_proc$run_id))
    if(NUM_RUNS_THIS_COMBO_ACTUAL_PROC == 0) next
    
    # --- Pre-cálculo general ---
    run_summary_info_proc <- raw_df_this_mean_sd_proc %>%
      group_by(run_id, social_distance_h, innovation_iul_Gamma) %>%
      summarise(
        adopters_prop_at_cell = num_adopters / N_nodes_actual,
        num_adopted_rational_at_cell = first(num_adopted_rational),
        num_adopted_social_at_cell = first(num_adopted_social),
        num_adopters_at_cell = first(num_adopters),
        initial_cluster_size_at_cell = first(initial_cluster_size),
        n_nodes_at_cell = first(N_nodes_actual),
        .groups = 'drop'
      ) %>%
      mutate(
        is_successful_at_cell = adopters_prop_at_cell >= SUCCESSFUL_DIFFUSION_THRESHOLD_PROP,
        new_adopters_at_cell = num_adopters_at_cell - initial_cluster_size_at_cell,
        prop_rational_of_new_at_cell = ifelse(new_adopters_at_cell > 0, num_adopted_rational_at_cell / new_adopters_at_cell, 0),
        prop_social_of_new_at_cell = ifelse(new_adopters_at_cell > 0, num_adopted_social_at_cell / new_adopters_at_cell, 0)
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
    
    # --- Llenar datos para los 4 heatmaps para esta (mean, sd) combo ---
    panel_data_tm_list_proc <- list()
    panel_data_aa_list_proc <- list()
    panel_data_pr_list_proc <- list()
    panel_data_ps_list_proc <- list()
    
    for (iul_val_proc in IUL_VALUES_SWEEP) {
      for (h_val_proc in H_VALUES_SWEEP) {
        runs_in_this_cell_proc <- run_summary_info_proc %>%
          filter(innovation_iul_Gamma == iul_val_proc, social_distance_h == h_val_proc)
        
        # Metric 1: Transition Metric
        prop_tm_cell <- NA_real_
        if (nrow(runs_in_this_cell_proc) > 0) {
          num_successful_cell_t <- sum(runs_in_this_cell_proc$is_successful_at_cell, na.rm = TRUE)
          if ((num_successful_cell_t / NUM_RUNS_THIS_COMBO_ACTUAL_PROC) >= MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL) {
            successful_and_transitioned_t <- runs_in_this_cell_proc %>%
              filter(is_successful_at_cell == TRUE & !is.na(first_transition_IUL_for_series) & first_transition_IUL_for_series <= iul_val_proc)
            prop_tm_cell <- nrow(successful_and_transitioned_t) / NUM_RUNS_THIS_COMBO_ACTUAL_PROC
          }
        }
        panel_data_tm_list_proc[[length(panel_data_tm_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_tm_cell)
        
        # Metric 2: Average Adoption
        prop_aa_cell <- NA_real_
        if(nrow(runs_in_this_cell_proc) > 0) {
          prop_aa_cell <- mean(runs_in_this_cell_proc$adopters_prop_at_cell, na.rm=TRUE)
        }
        panel_data_aa_list_proc[[length(panel_data_aa_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_aa_cell)
        
        # Metric 3: Proportion Rational
        prop_pr_cell <- NA_real_
        if(nrow(runs_in_this_cell_proc) > 0) {
          prop_pr_cell <- mean(runs_in_this_cell_proc$prop_rational_of_new_at_cell, na.rm=TRUE)
        }
        panel_data_pr_list_proc[[length(panel_data_pr_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_pr_cell)
        
        # Metric 4: Proportion Social
        prop_ps_cell <- NA_real_
        if(nrow(runs_in_this_cell_proc) > 0) {
          prop_ps_cell <- mean(runs_in_this_cell_proc$prop_social_of_new_at_cell, na.rm=TRUE)
        }
        panel_data_ps_list_proc[[length(panel_data_ps_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_ps_cell)
      }
    }
    # Add mean_tau parameter and combine
    df_tm_panel <- bind_rows(panel_data_tm_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, proportion_value_to_plot = val) %>% select(-val)
    df_aa_panel <- bind_rows(panel_data_aa_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, mean_adopters_prop_to_plot = val) %>% select(-val)
    df_pr_panel <- bind_rows(panel_data_pr_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, mean_prop_rational_to_plot = val) %>% select(-val)
    df_ps_panel <- bind_rows(panel_data_ps_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, mean_prop_social_to_plot = val) %>% select(-val)
    
    heatmap_data_tm_this_sd_list_proc[[mean_label_proc]] <- df_tm_panel
    heatmap_data_aa_this_sd_list_proc[[mean_label_proc]] <- df_aa_panel
    heatmap_data_pr_this_sd_list_proc[[mean_label_proc]] <- df_pr_panel
    heatmap_data_ps_this_sd_list_proc[[mean_label_proc]] <- df_ps_panel
    
  } # Fin bucle Medias de Tau para pre-procesamiento
  
  # Combinar todos los datos de las medias para el SD actual
  all_sds_transition_metric_heatmap_df_list[[sd_label_proc]] <- bind_rows(heatmap_data_tm_this_sd_list_proc[!sapply(heatmap_data_tm_this_sd_list_proc, is.null)])
  all_sds_avg_adoption_heatmap_df_list[[sd_label_proc]]    <- bind_rows(heatmap_data_aa_this_sd_list_proc[!sapply(heatmap_data_aa_this_sd_list_proc, is.null)])
  all_sds_prop_rational_heatmap_df_list[[sd_label_proc]]   <- bind_rows(heatmap_data_pr_this_sd_list_proc[!sapply(heatmap_data_pr_this_sd_list_proc, is.null)])
  all_sds_prop_social_heatmap_df_list[[sd_label_proc]]     <- bind_rows(heatmap_data_ps_this_sd_list_proc[!sapply(heatmap_data_ps_this_sd_list_proc, is.null)])
  
  cat(paste0("  Pre-procesamiento de datos para SD = ", current_tau_sd_proc, " completado.\n"))
}
cat("\nPre-procesamiento de todos los datos para heatmaps completado.\n")


# --- Función Auxiliar para crear UN heatmap individual (ya la tenías) ---
create_single_heatmap <- function(df_plot_data, fill_variable_name, legend_title_short, color_option="plasma", 
                                  show_legend=TRUE, y_axis_label_on=TRUE, plot_panel_title="") {
  if (is.null(df_plot_data) || nrow(df_plot_data) == 0 || all(is.na(df_plot_data[[fill_variable_name]]))) { # Modificado para chequear all(is.na())
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No plottable data\nfor this panel") + theme_void() + labs(title = plot_panel_title))
  }
  
  # Renombrar columnas iul y h si vienen con esos nombres del pre-procesamiento
  if("iul" %in% names(df_plot_data)) df_plot_data <- rename(df_plot_data, innovation_iul_Gamma = iul)
  if("h" %in% names(df_plot_data)) df_plot_data <- rename(df_plot_data, social_distance_h = h)
  
  df_plot_data$social_distance_h_factor <- factor(sprintf("%.2f", df_plot_data$social_distance_h), levels = sprintf("%.2f", H_VALUES_SWEEP))
  
  plot_obj <- ggplot(df_plot_data, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = .data[[fill_variable_name]])) +
    geom_tile(color = "grey90", lwd = 0.05) +
    scale_fill_viridis_c(name = legend_title_short, limits = c(0, 1), option=color_option, n.breaks=4, na.value = "white") +
    labs(x = NULL, y = if(y_axis_label_on) "MSP (h)" else NULL, title = plot_panel_title) +
    scale_x_continuous(breaks = seq(0, 1, 0.5), labels=c("0","0.5","1"), expand = c(0,0)) +
    scale_y_discrete(drop = FALSE) +
    theme_minimal(base_size = 7) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size=8),
      axis.text.x = element_text(angle = 45, hjust = 1, size=6), 
      axis.text.y = element_text(size=6), 
      axis.title.y = element_text(size=7, margin = margin(r = 0)),
      legend.position = if(show_legend) "right" else "none",
      legend.title = element_text(size = 7), legend.text = element_text(size = 6),
      legend.key.size = unit(0.6, "lines"),
      panel.grid = element_line(linewidth = 0.1)
    )
  return(plot_obj)
}


# --- Bucle Principal de PLOTEO: Generar un PDF por cada SD_tau ---
for (current_tau_sd_plot in TAU_NORMAL_SD_SWEEP_LIST) {
  sd_label_plot <- paste0("sd_", sprintf("%.2f", current_tau_sd_plot))
  cat(paste0("\nGenerating Consolidated PDF for Tau SD = ", current_tau_sd_plot, "\n"))
  
  # Extraer los dataframes de heatmap (ya procesados) para el SD actual
  df_tm_this_sd_plot <- all_sds_transition_metric_heatmap_df_list[[sd_label_plot]]
  df_aa_this_sd_plot <- all_sds_avg_adoption_heatmap_df_list[[sd_label_plot]]
  df_pr_this_sd_plot <- all_sds_prop_rational_heatmap_df_list[[sd_label_plot]]
  df_ps_this_sd_plot <- all_sds_prop_social_heatmap_df_list[[sd_label_plot]]
  
  if (is.null(df_tm_this_sd_plot) || is.null(df_aa_this_sd_plot) || is.null(df_pr_this_sd_plot) || is.null(df_ps_this_sd_plot)) {
    cat(paste0("  Faltan datos de heatmap procesados para SD = ", current_tau_sd_plot, ". Saltando PDF.\n"))
    next
  }
  
  plot_list_for_this_sd_pdf_final <- list()
  plot_titles_ordered <- c("Trans. Metric", "Avg Adoption", "Prop. Rational", "Prop. Social")
  fill_vars_ordered <- c("proportion_value_to_plot", "mean_adopters_prop_to_plot", "mean_prop_rational_to_plot", "mean_prop_social_to_plot")
  legend_titles_ordered <- c("Prop. Runs:\nSuccess & Trans.", "Mean Final\nAdopt. Prop.", "Mean Prop. New:\nRational", "Mean Prop. New:\nSocial Infl.")
  color_options_ordered <- c("plasma", "cividis", "magma", "inferno")
  data_frames_ordered <- list(df_tm_this_sd_plot, df_aa_this_sd_plot, df_pr_this_sd_plot, df_ps_this_sd_plot)
  
  for (current_threshold_mean_plot in THRESHOLD_MEAN_SWEEP_LIST) {
    for (j_metric in 1:4) { # Iterar sobre los 4 tipos de métricas/heatmaps
      
      current_df_for_plot <- data_frames_ordered[[j_metric]] %>% 
        filter(tau_mean_param == current_threshold_mean_plot)
      
      panel_title <- paste0("Mean τ=", sprintf("%.2f", current_threshold_mean_plot), "\n", plot_titles_ordered[j_metric])
      
      # Mostrar etiqueta Y solo para la primera columna de heatmaps (j_metric == 1)
      y_label_on <- (j_metric == 1)
      # Mostrar leyenda solo para la última columna de heatmaps (j_metric == 4)
      legend_on <- (j_metric == 4)
      
      plot_list_for_this_sd_pdf_final[[ ( (match(current_threshold_mean_plot, THRESHOLD_MEAN_SWEEP_LIST)-1)*4 ) + j_metric ]] <- 
        create_single_heatmap(
          current_df_for_plot, 
          fill_variable_name = fill_vars_ordered[j_metric], 
          legend_title_short = legend_titles_ordered[j_metric], 
          color_option = color_options_ordered[j_metric],
          show_legend = legend_on,
          y_axis_label_on = y_label_on,
          plot_panel_title = panel_title
        )
    }
  } 
  
  if (length(plot_list_for_this_sd_pdf_final) == 16) {
    combined_plot_for_pdf_final <- wrap_plots(plot_list_for_this_sd_pdf_final, ncol = 4, byrow = TRUE)
    
    # NUM_RUNS_SUBTITLE será el de la última combinación procesada, debería ser consistente
    num_runs_subtitle_val <- unique(raw_df_this_mean_sd_proc$run_id) %>% length()
    
    
    final_plot_with_annotation <- combined_plot_for_pdf_final + 
      plot_annotation(
        title = paste("Consolidated Heatmaps for ATP-net - Tau SD =", sprintf("%.2f", current_tau_sd_plot)),
        subtitle = paste("Thresholds ~ N(μ=var, σ=", sprintf("%.2f", current_tau_sd_plot), "). TM filter: ", MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL*100,"% runs in cell successful. ", num_runs_subtitle_val, " runs per (Γ,h) per individual panel.", sep=""),
        theme = theme(plot.title = element_text(hjust = 0.5, face="bold", size=14),
                      plot.subtitle = element_text(hjust = 0.5, size=9))
      ) 
    
    pdf_width <- 15 
    pdf_height <- 11 
    
    plot_filename_consolidated_final <- paste0(PLOTS_DIR, "consolidated_ALL_METRICS_sd", sprintf("%.2f", current_tau_sd_plot), ".pdf")
    ggsave(plot_filename_consolidated_final, final_plot_with_annotation, width = pdf_width, height = pdf_height, limitsize = FALSE)
    cat(paste0("  Saved consolidated PDF (All Metrics): ", plot_filename_consolidated_final, "\n"))
    
  } else {
    cat(paste0("  No se generaron suficientes plots (", length(plot_list_for_this_sd_pdf_final), "/16) para el PDF consolidado de SD = ", current_tau_sd_plot, ".\n"))
  }
} 

cat("\nGeneración de todos los PDFs consolidados (All Metrics) completada.\n")