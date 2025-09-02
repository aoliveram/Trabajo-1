# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script para Analizar Resultados CRUDOS y Generar PDF Consolidados con 16 Heatmaps
# Heatmap 2 (Phase Trans Prob) NO ACUMULATIVO, sensible a Γ y h.
# CORREGIDO: Transiciones en Γ=0 y h=0 no se cuentan.
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
H_VALUES_SWEEP <- seq(0/12, 12/12, by = 1/12) # Son 13 valores
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6)
TAU_NORMAL_SD_SWEEP_LIST <- c(0.08, 0.12, 0.16, 0.20)

PHASE_TRANSITION_THRESHOLD_JUMP <- 1/3 # NUEVO VALOR
# MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL ya no se usa para el heatmap de transiciones
# SUCCESSFUL_DIFFUSION_THRESHOLD_PROP (0.50) se usa en el pre-cálculo general para la columna is_successful_at_cell

cat("Cargando resultados crudos guardados...\n")
# grand_raw_results_path <- paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_raw_results_all_sds_means_2.rds") # Nombre anterior
grand_raw_results_path <- paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_raw_results_all_sds_means.rds") # NUEVO NOMBRE
if (!file.exists(grand_raw_results_path)) {
  stop("Archivo de resultados crudos no encontrado: ", grand_raw_results_path)
}
all_sds_raw_results_list_from_file <- readRDS(grand_raw_results_path) # Renombrar para evitar confusión con el del script anterior
cat("Resultados cargados.\n")

# --- PASO 0: Pre-procesar los datos crudos ---
all_sds_transition_metric_heatmap_df_list <- list()
all_sds_avg_adoption_heatmap_df_list    <- list()
all_sds_prop_rational_heatmap_df_list   <- list()
all_sds_prop_social_heatmap_df_list     <- list()

cat("\nPre-procesando datos crudos para todos los heatmaps...\n")
for (current_tau_sd_proc in TAU_NORMAL_SD_SWEEP_LIST) {
  sd_label_proc <- paste0("sd_", sprintf("%.2f", current_tau_sd_proc))
  current_sd_all_means_raw_results_proc <- all_sds_raw_results_list_from_file[[sd_label_proc]] # Usar el objeto cargado
  if (is.null(current_sd_all_means_raw_results_proc)) { next }
  
  heatmap_data_tm_this_sd_list_proc <- list(); heatmap_data_aa_this_sd_list_proc <- list()
  heatmap_data_pr_this_sd_list_proc <- list(); heatmap_data_ps_this_sd_list_proc <- list()
  
  for (current_threshold_mean_proc in THRESHOLD_MEAN_SWEEP_LIST) {
    mean_label_proc <- paste0("mean_", sprintf("%.2f", current_threshold_mean_proc))
    raw_df_this_mean_sd_proc <- current_sd_all_means_raw_results_proc[[mean_label_proc]]
    if (is.null(raw_df_this_mean_sd_proc) || nrow(raw_df_this_mean_sd_proc) == 0) { next }
    
    NUM_RUNS_THIS_COMBO_ACTUAL_PROC <- length(unique(raw_df_this_mean_sd_proc$run_id))
    if(NUM_RUNS_THIS_COMBO_ACTUAL_PROC == 0) next
    
    base_run_summary_proc <- raw_df_this_mean_sd_proc %>%
      group_by(run_id, social_distance_h, innovation_iul_Gamma) %>% 
      summarise(
        adopters_prop_at_cell = first(num_adopters) / first(N_nodes_actual),
        num_adopted_rational_at_cell = first(num_adopted_rational),
        num_adopted_social_at_cell = first(num_adopted_social),
        num_adopters_at_cell = first(num_adopters),
        initial_cluster_size_at_cell = first(initial_cluster_size),
        .groups = 'drop'
      ) %>%
      mutate(
        new_adopters_at_cell = num_adopters_at_cell - initial_cluster_size_at_cell,
        prop_rational_of_new_at_cell = ifelse(new_adopters_at_cell > 0, num_adopted_rational_at_cell / new_adopters_at_cell, 0),
        prop_social_of_new_at_cell = ifelse(new_adopters_at_cell > 0, num_adopted_social_at_cell / new_adopters_at_cell, 0)
      )
    
    # 1. Identificar primeras transiciones por CAMBIO EN GAMMA
    first_transitions_by_gamma <- base_run_summary_proc %>%
      arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
      group_by(run_id, social_distance_h) %>%
      mutate(
        # El lag default es NA. Si es el primer Γ, el salto será NA.
        jump_vs_prev_gamma = adopters_prop_at_cell - lag(adopters_prop_at_cell), 
        is_transition_vs_prev_gamma = ifelse(!is.na(jump_vs_prev_gamma) & jump_vs_prev_gamma >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
      ) %>%
      filter(is_transition_vs_prev_gamma == 1) %>%
      summarise(first_transition_IUL = min(innovation_iul_Gamma), .groups = 'drop')
    
    # 2. Identificar primeras transiciones por CAMBIO EN H
    first_transitions_by_h <- base_run_summary_proc %>%
      arrange(run_id, innovation_iul_Gamma, social_distance_h) %>% 
      group_by(run_id, innovation_iul_Gamma) %>%
      mutate(
        # El lag default es NA. Si es el primer h, el salto será NA.
        jump_vs_prev_h = adopters_prop_at_cell - lag(adopters_prop_at_cell),
        is_transition_vs_prev_h = ifelse(!is.na(jump_vs_prev_h) & jump_vs_prev_h >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
      ) %>%
      filter(is_transition_vs_prev_h == 1) %>%
      summarise(first_transition_H = min(social_distance_h), .groups = 'drop')
    
    panel_data_tm_list_proc <- list(); panel_data_aa_list_proc <- list()
    panel_data_pr_list_proc <- list(); panel_data_ps_list_proc <- list()
    
    for (iul_val_proc in IUL_VALUES_SWEEP) {
      for (h_val_proc in H_VALUES_SWEEP) {
        
        # Metric 1: Transition Metric (NO ACUMULATIVA, sensible a Γ y h, SIN filtro de éxito de celda)
        runs_transitioned_to_this_cell_by_gamma <- first_transitions_by_gamma %>%
          filter(social_distance_h == h_val_proc, first_transition_IUL == iul_val_proc) %>%
          pull(run_id)
        runs_transitioned_to_this_cell_by_h <- first_transitions_by_h %>%
          filter(innovation_iul_Gamma == iul_val_proc, first_transition_H == h_val_proc) %>%
          pull(run_id)
        unique_runs_transitioned_to_cell <- unique(c(runs_transitioned_to_this_cell_by_gamma, runs_transitioned_to_this_cell_by_h))
        prop_tm_cell_non_cumulative_combined <- length(unique_runs_transitioned_to_cell) / NUM_RUNS_THIS_COMBO_ACTUAL_PROC
        
        panel_data_tm_list_proc[[length(panel_data_tm_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_tm_cell_non_cumulative_combined)
        
        # Métricas 2, 3, 4
        current_cell_data_proc <- base_run_summary_proc %>%
          filter(innovation_iul_Gamma == iul_val_proc, social_distance_h == h_val_proc)
        prop_aa_cell <- NA_real_; if(nrow(current_cell_data_proc) == NUM_RUNS_THIS_COMBO_ACTUAL_PROC) {prop_aa_cell <- mean(current_cell_data_proc$adopters_prop_at_cell, na.rm=TRUE); if(is.nan(prop_aa_cell)) prop_aa_cell <- NA_real_}
        panel_data_aa_list_proc[[length(panel_data_aa_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_aa_cell)
        prop_pr_cell <- NA_real_; if(nrow(current_cell_data_proc) == NUM_RUNS_THIS_COMBO_ACTUAL_PROC) {prop_pr_cell <- mean(current_cell_data_proc$prop_rational_of_new_at_cell, na.rm=TRUE); if(is.nan(prop_pr_cell)) prop_pr_cell <- NA_real_}
        panel_data_pr_list_proc[[length(panel_data_pr_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_pr_cell)
        prop_ps_cell <- NA_real_; if(nrow(current_cell_data_proc) == NUM_RUNS_THIS_COMBO_ACTUAL_PROC) {prop_ps_cell <- mean(current_cell_data_proc$prop_social_of_new_at_cell, na.rm=TRUE); if(is.nan(prop_ps_cell)) prop_ps_cell <- NA_real_}
        panel_data_ps_list_proc[[length(panel_data_ps_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_ps_cell)
      }
    }
    heatmap_data_tm_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_tm_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, proportion_value_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)
    heatmap_data_aa_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_aa_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, mean_adopters_prop_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)
    heatmap_data_pr_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_pr_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, mean_prop_rational_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)
    heatmap_data_ps_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_ps_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, mean_prop_social_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)
  }
  all_sds_transition_metric_heatmap_df_list[[sd_label_proc]] <- bind_rows(heatmap_data_tm_this_sd_list_proc[!sapply(heatmap_data_tm_this_sd_list_proc, is.null)])
  all_sds_avg_adoption_heatmap_df_list[[sd_label_proc]]    <- bind_rows(heatmap_data_aa_this_sd_list_proc[!sapply(heatmap_data_aa_this_sd_list_proc, is.null)])
  all_sds_prop_rational_heatmap_df_list[[sd_label_proc]]   <- bind_rows(heatmap_data_pr_this_sd_list_proc[!sapply(heatmap_data_pr_this_sd_list_proc, is.null)])
  all_sds_prop_social_heatmap_df_list[[sd_label_proc]]     <- bind_rows(heatmap_data_ps_this_sd_list_proc[!sapply(heatmap_data_ps_this_sd_list_proc, is.null)])
}
cat("Pre-procesamiento de todos los datos para heatmaps completado (con TM no acumulativa y sensible a Γ/h).\n")

# --- Función Auxiliar para crear UN heatmap individual ---
# (Misma función `create_single_heatmap_v2` que en la respuesta anterior)
create_single_heatmap_v2 <- function(df_plot_data, fill_col_name, legend_title, viridis_option, 
                                     show_legend=TRUE, y_axis_label_on=TRUE, x_axis_label_on=TRUE, 
                                     panel_super_title="") {
  if (is.null(df_plot_data) || nrow(df_plot_data) == 0 || all(is.na(df_plot_data[[fill_col_name]]))) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No plottable data") + theme_void() + labs(title = panel_super_title))
  }
  if("iul" %in% names(df_plot_data)) df_plot_data <- rename(df_plot_data, innovation_iul_Gamma = iul)
  if("h" %in% names(df_plot_data)) df_plot_data <- rename(df_plot_data, social_distance_h = h)
  df_plot_data$social_distance_h_factor <- factor(sprintf("%.2f", df_plot_data$social_distance_h), levels = sprintf("%.2f", H_VALUES_SWEEP))
  plot_obj <- ggplot(df_plot_data, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = .data[[fill_col_name]])) +
    geom_tile(color = "grey90", lwd = 0.05) +
    scale_fill_viridis_c(name = legend_title, limits = c(0, 1), option=viridis_option, n.breaks=4, na.value = "white") +
    labs(x = if(x_axis_label_on) expression(paste("IUL (", Gamma, ")")) else NULL, y = if(y_axis_label_on) "MSP (h)" else NULL, title = panel_super_title) +
    scale_x_continuous(breaks = seq(0, 1, 0.5), labels= if(x_axis_label_on) c("0","0.5","1") else waiver(), expand = c(0,0)) +
    scale_y_discrete(drop = FALSE, labels = if(y_axis_label_on) waiver() else NULL ) + 
    theme_minimal(base_size = 6) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=7),
          axis.text.x = element_text(angle = 45, hjust = 1, size=5), axis.text.y = element_text(size=5), 
          axis.title.x = element_text(size=6, margin = margin(t = 1)), 
          axis.title.y = element_text(size=6, margin = margin(r = 1)), 
          legend.position = if(show_legend) "right" else "none",
          legend.title = element_text(size = 6), legend.text = element_text(size = 5),
          legend.key.size = unit(0.5, "lines"), panel.grid = element_line(linewidth = 0.1),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")) 
  return(plot_obj)
}

# --- Bucle Principal de PLOTEO: Generar un PDF por cada MEAN_tau ---
# (Misma lógica de ploteo que en la respuesta anterior, 
#  solo el título de la leyenda para el heatmap 2 y el subtítulo general cambiarán)
for (current_threshold_mean_plot in THRESHOLD_MEAN_SWEEP_LIST) {
  cat(paste0("\nGenerating Consolidated PDF for Mean τ = ", current_threshold_mean_plot, "\n"))
  plot_list_for_this_mean_pdf_ordered <- vector("list", length(TAU_NORMAL_SD_SWEEP_LIST) * 4)
  
  metric_titles_for_cols <- c("Avg. Adoption", "Phase Trans Prob.\n(Exact Γ or h Jump)", "Prop. Rational", "Prop. Social Infl.")
  metric_fill_vars_for_cols <- c("mean_adopters_prop_to_plot", "proportion_value_to_plot", "mean_prop_rational_to_plot", "mean_prop_social_to_plot")
  metric_legend_titles_for_cols <- c("Mean Final\nAdopt. Prop.", "Prop. Runs with\n1st Transition\nEXACTLY at this\n(Γ or h) Point", # LEYENDA CAMBIADA
                                     "Mean Prop. New:\nRational", "Mean Prop. New:\nSocial Infl.")
  metric_color_options_for_cols <- rep("viridis", 4)
  
  data_sources_for_cols <- list(
    all_sds_avg_adoption_heatmap_df_list,
    all_sds_transition_metric_heatmap_df_list, # Contiene la nueva métrica
    all_sds_prop_rational_heatmap_df_list,
    all_sds_prop_social_heatmap_df_list
  )
  
  for (row_idx in 1:length(TAU_NORMAL_SD_SWEEP_LIST)) {
    current_tau_sd_plot <- TAU_NORMAL_SD_SWEEP_LIST[row_idx]
    sd_label_plot <- paste0("sd_", sprintf("%.2f", current_tau_sd_plot))
    for (col_idx in 1:4) { 
      df_all_means_for_sd_metric <- data_sources_for_cols[[col_idx]][[sd_label_plot]]
      current_df_for_panel <- if(!is.null(df_all_means_for_sd_metric)) {
        df_all_means_for_sd_metric %>% filter(tau_mean_param == current_threshold_mean_plot)
      } else { data.frame() }
      panel_individual_title <- if(col_idx == 1) paste0("SD τ=",sprintf("%.2f", current_tau_sd_plot)) else ""
      y_label_on <- (col_idx == 1); x_label_on <- (row_idx == length(TAU_NORMAL_SD_SWEEP_LIST)); legend_on <- (col_idx == 4)
      plot_index_in_list <- ( (row_idx-1)*4 ) + col_idx
      plot_list_for_this_mean_pdf_ordered[[plot_index_in_list]] <- create_single_heatmap_v2(
        current_df_for_panel, metric_fill_vars_for_cols[col_idx], metric_legend_titles_for_cols[col_idx], 
        metric_color_options_for_cols[col_idx], legend_on, y_label_on, x_label_on, panel_individual_title)
    }
  } 
  
  if (length(plot_list_for_this_mean_pdf_ordered) == (length(TAU_NORMAL_SD_SWEEP_LIST) * 4) ) {
    col_titles_plots <- lapply(metric_titles_for_cols, function(title) {ggplot() + labs(title=title) + theme_void() + theme(plot.title = element_text(hjust=0.5, size=9, face="bold", margin = margin(b=1, unit="mm")))})
    column_titles_row_layout <- Reduce(`+`, col_titles_plots) + plot_layout(ncol = 4)
    heatmaps_grid_layout <- wrap_plots(plot_list_for_this_mean_pdf_ordered, ncol = 4, byrow = TRUE)
    final_combined_layout <- column_titles_row_layout / heatmaps_grid_layout + plot_layout(heights = c(0.07, 1)) 
    
    num_runs_val_subtitle <- NA 
    first_valid_sd_label <- names(all_sds_raw_results_list_from_file)[which(!sapply(all_sds_raw_results_list_from_file, is.null))[1]]
    if(!is.na(first_valid_sd_label) && !is.null(all_sds_raw_results_list_from_file[[first_valid_sd_label]])){
      first_valid_mean_label <- names(all_sds_raw_results_list_from_file[[first_valid_sd_label]])[which(!sapply(all_sds_raw_results_list_from_file[[first_valid_sd_label]], is.null))[1]]
      if(!is.na(first_valid_mean_label) && !is.null(all_sds_raw_results_list_from_file[[first_valid_sd_label]][[first_valid_mean_label]])){
        num_runs_val_subtitle <- length(unique(all_sds_raw_results_list_from_file[[first_valid_sd_label]][[first_valid_mean_label]]$run_id))
      }
    }
    if(is.na(num_runs_val_subtitle)) num_runs_val_subtitle <- "N/A"
    
    
    final_plot_with_annotation <- final_combined_layout + 
      plot_annotation(
        title = paste("Consolidated Heatmaps for ATP-net - Mean τ =", sprintf("%.2f", current_threshold_mean_plot)),
        # Subtítulo ahora no menciona el filtro de TM ya que no aplica a todos los heatmaps de la misma manera
        subtitle = paste("Thresholds ~ N(μ=", sprintf("%.2f", current_threshold_mean_plot), ", σ=var). ", num_runs_val_subtitle, " runs per (Γ,h) per individual panel.", sep=""),
        theme = theme(plot.title = element_text(hjust = 0.5, face="bold", size=14), plot.subtitle = element_text(hjust = 0.5, size=8))
      ) 
    pdf_width <- 12; pdf_height <- 13 
    plot_filename_consolidated_final <- paste0(PLOTS_DIR, "consolidated_ALL_METRICS_TM_CORRECTED_mean_tau_", sprintf("%.2f", current_threshold_mean_plot), ".pdf")
    ggsave(plot_filename_consolidated_final, final_plot_with_annotation, width = pdf_width, height = pdf_height, limitsize = FALSE)
    cat(paste0("  Saved consolidated PDF (Corrected TM): ", plot_filename_consolidated_final, "\n"))
  } else {
    cat(paste0("  No se generaron suficientes plots para el PDF consolidado de Mean τ = ", current_threshold_mean_plot, ".\n"))
  }
} 
cat("\nGeneración de todos los PDFs consolidados (Corrected TM) completada.\n")

saveRDS(all_sds_transition_metric_heatmap_df_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_CORRECTED_TM_heatmap_data_all_sds.rds"))
cat("Datos para heatmaps (Corrected TM) guardados.\n")