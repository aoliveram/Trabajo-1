# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script para Analizar Resultados raw y generar PDF con 16 Heatmaps
#
# Entradas (RDS crudos):
#   - trabajo_1_files/GSS_diffusion_simulation_files_sigm/GSS_phase_transition_GRAND_COMBINED... .rds
#
# Umbrales/convenciones:
#   - PHASE_TRANSITION_THRESHOLD_JUMP = 1/3 (salto mínimo para marcar transición).
#   - SUCCESSFUL_DIFFUSION_THRESHOLD_PROP = 0.50 (no se usa como filtro de celda para TM en este script).
#
# Plots:
#   1) Avg. Adoption (proporción media TOTAL de adoptantes por celda).
#   2) Phase Trans Prob. (No acumulativa y sin filtro de "celda exitosa").
#   3) Avg. Adopt. by Rational (of population): media de num_adopted_rational / N_nodes_actual.
#   4) Avg. Adopt. by Social (of population): media de num_adopted_social / N_nodes_actual.
#
#   - Nombre de salida: PLOTS_DIR/heatmaps__seed_<SEEDING_STRATEGY_FIXED>_mean_tau_<XX>.pdf
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- 0. Cargar Librerías ---
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(viridis)

# --- 1. Parámetros de Análisis ---

# ---- MODIFICACIÓN 1: Cambiar etiquetas y apuntar a los archivos de GSS ----
RESULTS_DIR <- "trabajo_1_files/GSS_diffusion_simulation_files_sigm/"
PLOTS_DIR <- "trabajo_1_plots/GSS_diffusion_simulation_plots_sigm/"
NETWORK_LABEL <- "GSS-net" # Etiqueta para los títulos de los gráficos
# --------------------------------------------------------------------------

dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Parámetros (deben ser idénticos a los del script de simulación)
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP <- seq(0/12, 12/12, by = 1/12)
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6)
TAU_NORMAL_SD_SWEEP_LIST <- c(0.08, 0.12, 0.16, 0.20)

# Estrategia de siembra que se va a plotear. ¡Asegúrate que coincida con los archivos que quieres analizar!
strategies <- c("random", "central", "marginal", "eigen", "closeness")
SEEDING_STRATEGY_FIXED <- strategies[1] # Cambiar según sea necesario

# Parámetros de umbral para el análisis
PHASE_TRANSITION_THRESHOLD_JUMP <- 1/3
SUCCESSFUL_DIFFUSION_THRESHOLD_PROP <- 0.50

# --- 2. Carga y Pre-procesamiento de Datos ---

cat("Cargando resultados crudos de la simulación GSS...\n")
# ---- MODIFICACIÓN 2: Cambiar el nombre del archivo de resultados de ATP a GSS ----
grand_raw_results_path <- paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_raw_results_all_sds_means_", SEEDING_STRATEGY_FIXED, ".rds")

# ---------------------------------------------------------------------------------

if (!file.exists(grand_raw_results_path)) {
  stop("Archivo de resultados crudos no encontrado: ", grand_raw_results_path)
}
all_sds_raw_results_list_from_file <- readRDS(grand_raw_results_path)
cat("Resultados cargados.\n")

# Pre-procesamiento para agregar las métricas de interés
# (Esta sección es idéntica en lógica a la del script de ATP)
all_sds_transition_metric_heatmap_df_list <- list()
all_sds_avg_adoption_heatmap_df_list    <- list()
all_sds_avg_rational_adopt_pop_heatmap_df_list <- list()
all_sds_avg_social_adopt_pop_heatmap_df_list   <- list()

cat("\nPre-procesando datos crudos para generar los heatmaps...\n")
for (current_tau_sd_proc in TAU_NORMAL_SD_SWEEP_LIST) {
  sd_label_proc <- paste0("sd_", sprintf("%.2f", current_tau_sd_proc))
  current_sd_all_means_raw_results_proc <- all_sds_raw_results_list_from_file[[sd_label_proc]]
  if (is.null(current_sd_all_means_raw_results_proc)) { next }
  
  heatmap_data_tm_this_sd_list_proc <- list(); heatmap_data_aa_this_sd_list_proc <- list()
  heatmap_data_ar_pop_this_sd_list_proc <- list(); heatmap_data_as_pop_this_sd_list_proc <- list()
  
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
        prop_rational_of_pop_at_cell = first(num_adopted_rational) / first(N_nodes_actual),
        prop_social_of_pop_at_cell = first(num_adopted_social) / first(N_nodes_actual),
        .groups = 'drop'
      ) %>%
      arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
      group_by(run_id, social_distance_h) %>%
      mutate(
        jump_at_step = adopters_prop_at_cell - lag(adopters_prop_at_cell), 
        is_transition_vs_prev_gamma = ifelse(!is.na(jump_at_step) & jump_at_step >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
      ) %>%
      mutate(first_transition_IUL_for_series = if (any(is_transition_vs_prev_gamma == 1, na.rm=TRUE)) {min(innovation_iul_Gamma[which(is_transition_vs_prev_gamma == 1)])} else {NA_real_}) %>%
      ungroup() %>%
      arrange(run_id, innovation_iul_Gamma, social_distance_h) %>%
      group_by(run_id, innovation_iul_Gamma) %>%
      mutate(
        jump_vs_prev_h = adopters_prop_at_cell - lag(adopters_prop_at_cell),
        is_transition_vs_prev_h = ifelse(!is.na(jump_vs_prev_h) & jump_vs_prev_h >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
      ) %>%
      mutate(first_transition_H_for_series = if (any(is_transition_vs_prev_h == 1, na.rm=TRUE)) {min(social_distance_h[which(is_transition_vs_prev_h == 1)])} else {NA_real_}) %>%
      ungroup()
    
    panel_data_tm_list_proc <- list(); panel_data_aa_list_proc <- list()
    panel_data_ar_pop_list_proc <- list(); panel_data_as_pop_list_proc <- list()
    
    for (iul_val_proc in IUL_VALUES_SWEEP) {
      for (h_val_proc in H_VALUES_SWEEP) {
        
        # Metric 1: Phase Transition Probability
        runs_transitioned_to_this_cell_by_gamma <- base_run_summary_proc %>%
          filter(social_distance_h == h_val_proc, !is.na(first_transition_IUL_for_series), first_transition_IUL_for_series == iul_val_proc) %>% pull(run_id)
        runs_transitioned_to_this_cell_by_h <- base_run_summary_proc %>%
          filter(innovation_iul_Gamma == iul_val_proc, !is.na(first_transition_H_for_series), first_transition_H_for_series == h_val_proc) %>% pull(run_id)
        unique_runs_transitioned_to_cell <- unique(c(runs_transitioned_to_this_cell_by_gamma, runs_transitioned_to_this_cell_by_h))
        prop_tm_cell_combined <- length(unique_runs_transitioned_to_cell) / NUM_RUNS_THIS_COMBO_ACTUAL_PROC
        panel_data_tm_list_proc[[length(panel_data_tm_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_tm_cell_combined)
        
        # Métricas 2, 3, 4
        current_cell_data_proc <- base_run_summary_proc %>%
          filter(innovation_iul_Gamma == iul_val_proc, social_distance_h == h_val_proc)
        
        prop_aa_cell <- mean(current_cell_data_proc$adopters_prop_at_cell, na.rm=TRUE)
        panel_data_aa_list_proc[[length(panel_data_aa_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_aa_cell)
        
        prop_ar_pop_cell <- mean(current_cell_data_proc$prop_rational_of_pop_at_cell, na.rm=TRUE)
        panel_data_ar_pop_list_proc[[length(panel_data_ar_pop_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_ar_pop_cell)
        
        prop_as_pop_cell <- mean(current_cell_data_proc$prop_social_of_pop_at_cell, na.rm=TRUE)
        panel_data_as_pop_list_proc[[length(panel_data_as_pop_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_as_pop_cell)
      }
    }
    
    heatmap_data_tm_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_tm_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, proportion_value_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)
    heatmap_data_aa_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_aa_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, mean_adopters_prop_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)
    heatmap_data_ar_pop_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_ar_pop_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, avg_rational_adopt_pop_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)
    heatmap_data_as_pop_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_as_pop_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, avg_social_adopt_pop_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)
  }
  
  all_sds_transition_metric_heatmap_df_list[[sd_label_proc]] <- bind_rows(heatmap_data_tm_this_sd_list_proc[!sapply(heatmap_data_tm_this_sd_list_proc, is.null)])
  all_sds_avg_adoption_heatmap_df_list[[sd_label_proc]]    <- bind_rows(heatmap_data_aa_this_sd_list_proc[!sapply(heatmap_data_aa_this_sd_list_proc, is.null)])
  all_sds_avg_rational_adopt_pop_heatmap_df_list[[sd_label_proc]] <- bind_rows(heatmap_data_ar_pop_this_sd_list_proc[!sapply(heatmap_data_ar_pop_this_sd_list_proc, is.null)])
  all_sds_avg_social_adopt_pop_heatmap_df_list[[sd_label_proc]]   <- bind_rows(heatmap_data_as_pop_this_sd_list_proc[!sapply(heatmap_data_as_pop_this_sd_list_proc, is.null)])
}
cat("Pre-procesamiento de todos los datos para heatmaps completado.\n")


# --- 3. Función Genérica para Crear un Heatmap ---
# (Esta función no necesita cambios, es genérica)
create_single_heatmap <- function(df_plot_data, fill_col_name, legend_title_text, viridis_option, 
                                  show_legend=TRUE, y_axis_label_on=TRUE, x_axis_label_on=TRUE, 
                                  panel_row_title="") {
  
  if (is.null(df_plot_data) || nrow(df_plot_data) == 0 || all(is.na(df_plot_data[[fill_col_name]]))) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data") + theme_void() + 
             labs(y = if(y_axis_label_on) panel_row_title else NULL) + 
             theme(axis.title.y = element_text(size=8, face="bold", angle=90, margin = margin(r=5))))
  }
  
  df_plot_data$social_distance_h_factor <- factor(sprintf("%.2f", df_plot_data$h), levels = sprintf("%.2f", sort(unique(H_VALUES_SWEEP))))
  y_breaks <- sprintf("%.2f", H_VALUES_SWEEP[seq(1, length(H_VALUES_SWEEP), by = 2)])
  x_breaks <- seq(0, 1, 0.25)
  
  ggplot(df_plot_data, aes(x = iul, y = social_distance_h_factor, fill = .data[[fill_col_name]])) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(name = if(show_legend) legend_title_text else NULL, limits = c(0, 1), option=viridis_option, n.breaks=4, na.value = "grey90") +
    labs(
      x = if(x_axis_label_on) expression(paste("IUL (", Gamma, ")")) else NULL, 
      y = if(y_axis_label_on) panel_row_title else NULL,
      title = NULL
    ) +
    scale_x_continuous(breaks = x_breaks, labels = if(x_axis_label_on) sprintf("%.2f", x_breaks) else NULL, expand = c(0,0)) +
    scale_y_discrete(drop = FALSE, breaks = y_breaks, labels = if(y_axis_label_on) y_breaks else NULL) + 
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=6, color = if(x_axis_label_on) "black" else "transparent"),
      axis.text.y = element_text(size=6, color = if(y_axis_label_on) "black" else "transparent"),
      axis.title.x = element_text(size=8, face="bold", margin = margin(t = 2, unit="mm")),
      axis.title.y = element_text(size=8, face="bold", angle=90, margin = margin(r = 2, unit="mm")),
      legend.position = if(show_legend) "right" else "none",
      legend.title = element_text(size = 7), 
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.6, "lines"), 
      panel.grid = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")
    )
}

# --- 4. Bucle Principal de Ploteo: Generar un PDF por cada Media de Umbral ---
for (current_threshold_mean_plot in THRESHOLD_MEAN_SWEEP_LIST) {
  cat(paste0("\nGenerando PDF consolidado para Media del Umbral μ = ", current_threshold_mean_plot, "\n"))
  
  plot_list_for_this_mean_pdf <- list()
  
  metric_titles_for_cols <- c("Avg. Adoption", "Phase Trans Prob.", "Avg. Adopt. by Rational", "Avg. Adopt. by Social Infl.")
  metric_fill_vars_for_cols <- c("mean_adopters_prop_to_plot", "proportion_value_to_plot", "avg_rational_adopt_pop_to_plot", "avg_social_adopt_pop_to_plot")
  data_sources_for_cols <- list(all_sds_avg_adoption_heatmap_df_list, all_sds_transition_metric_heatmap_df_list, all_sds_avg_rational_adopt_pop_heatmap_df_list, all_sds_avg_social_adopt_pop_heatmap_df_list)
  
  for (row_idx in 1:length(TAU_NORMAL_SD_SWEEP_LIST)) {
    current_tau_sd_plot <- TAU_NORMAL_SD_SWEEP_LIST[row_idx]
    sd_label_plot <- paste0("sd_", sprintf("%.2f", current_tau_sd_plot))
    
    for (col_idx in 1:4) {
      df_all_means_for_sd_metric <- data_sources_for_cols[[col_idx]][[sd_label_plot]]
      
      current_df_for_panel <- if (!is.null(df_all_means_for_sd_metric) && nrow(df_all_means_for_sd_metric) > 0) {
        df_all_means_for_sd_metric %>% filter(tau_mean_param == current_threshold_mean_plot)
      } else {
        data.frame()
      }
      
      row_title_str <- if (col_idx == 1) paste0("MSP (h) - SD=", sprintf("%.2f", current_tau_sd_plot)) else ""
      y_label_visible <- (col_idx == 1)
      x_label_visible <- (row_idx == length(TAU_NORMAL_SD_SWEEP_LIST))
      legend_visible <- (col_idx == 4)
      
      plot_list_for_this_mean_pdf[[ ( (row_idx-1)*4 ) + col_idx ]] <- 
        create_single_heatmap(
          df_plot_data = current_df_for_panel, 
          fill_col_name = metric_fill_vars_for_cols[col_idx], 
          legend_title_text = NULL,
          viridis_option = "viridis",
          show_legend = legend_visible,
          y_axis_label_on = y_label_visible,
          x_axis_label_on = x_label_visible,
          panel_row_title = row_title_str
        )
    }
  }
  
  if (length(plot_list_for_this_mean_pdf) == (length(TAU_NORMAL_SD_SWEEP_LIST) * 4)) {
    
    col_titles_plots <- lapply(metric_titles_for_cols, function(title) {
      ggplot() + labs(title=title) + theme_void() + 
        theme(plot.title = element_text(hjust=0.5, size=10, face="bold", margin = margin(b=0, t=2, unit="mm")))
    })
    
    column_titles_row_layout <- wrap_plots(col_titles_plots, ncol = 4)
    heatmaps_grid_layout <- wrap_plots(plot_list_for_this_mean_pdf, ncol = 4, byrow = TRUE)
    
    final_combined_layout <- column_titles_row_layout / heatmaps_grid_layout + 
      plot_layout(heights = c(0.05, 1))
    
    num_runs_val_subtitle <- NA
    try({
      num_runs_val_subtitle <- length(unique(data_sources_for_cols[[1]][[1]][[1]]$run_id))
    }, silent = TRUE)
    if(is.na(num_runs_val_subtitle)) num_runs_val_subtitle <- "N/A"
    
    # ---- MODIFICACIÓN 6: Cambiar el título del PDF de ATP a GSS ----
    final_plot_with_annotation <- final_combined_layout + 
      plot_annotation(
        title = paste("Consolidated Heatmaps for", NETWORK_LABEL, "- Mean threshold =", sprintf("%.2f", current_threshold_mean_plot)),
        subtitle = paste("Thresholds ~ N(μ=", sprintf("%.2f", current_threshold_mean_plot), ", SD=var). ", 
                         num_runs_val_subtitle, " runs per cell. Seeding strategy: ", SEEDING_STRATEGY_FIXED,
                         sep=""),
        theme = theme(plot.title = element_text(hjust = 0.5, face="bold", size=12),
                      plot.subtitle = element_text(hjust = 0.5, size=9.5))
      )
    # -----------------------------------------------------------------
    
    # ---- MODIFICACIÓN 7: Cambiar el nombre del archivo PDF de ATP a GSS ----
    plot_filename <- file.path(PLOTS_DIR, paste0("heatmaps_GSS_seed_", SEEDING_STRATEGY_FIXED, "_mean_tau_", sprintf("%.2f", current_threshold_mean_plot), ".pdf"))
    ggsave(plot_filename, final_plot_with_annotation, width = 7.5, height = 7.0, limitsize = FALSE)
    cat(paste0("  PDF consolidado guardado en: ", plot_filename, "\n"))
    # ----------------------------------------------------------------------
    
  } else {
    cat(paste0("  No se generaron suficientes plots para el PDF de μ = ", current_threshold_mean_plot, ".\n"))
  }
} 
cat("\nGeneración de todos los PDFs finalizada.\n")