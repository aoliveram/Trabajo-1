# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script para Analizar Resultados raw y generar PDF con 16 Heatmaps (SD)
#
# Salida:
#   - 4 columnas: SD(Adopción Total), SD(Transición), SD(Adopción Racional), SD(Adopción Social)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(viridis)

# --- Parámetros de Análisis ---
RESULTS_DIR <- "trabajo_1_files/ATP_diffusion_simulation_files_sigm/"
PLOTS_DIR <- "trabajo_1_plots/ATP_diffusion_simulation_plots_sigm/"
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP <- seq(0/12, 12/12, by = 1/12)
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6)
TAU_NORMAL_SD_SWEEP_LIST <- c(0.08, 0.12, 0.16, 0.20)

strategies <- c("random", "central", "marginal", "eigen", "closeness")
SEEDING_STRATEGY_FIXED <- strategies[2]
PHASE_TRANSITION_THRESHOLD_JUMP <- 1/3

cat("Cargando resultados crudos guardados...\n")
grand_raw_results_path <- paste0(
  RESULTS_DIR,
  "phase_transition_GRAND_COMBINED_raw_results_all_sds_means_",
  SEEDING_STRATEGY_FIXED, ".rds"
)
if (!file.exists(grand_raw_results_path)) {
  stop("Archivo de resultados crudos no encontrado: ", grand_raw_results_path)
}
all_sds_raw_results_list_from_file <- readRDS(grand_raw_results_path)
cat("Resultados cargados.\n")

# --- Inicializar valor máximo global de SD ---
max_sd_global <- 0

# --- PASO 0: Pre-procesar los datos crudos (para SD) ---
all_sds_transition_metric_heatmap_df_list <- list()
all_sds_avg_adoption_heatmap_df_list <- list()
all_sds_avg_rational_adopt_pop_heatmap_df_list <- list()
all_sds_avg_social_adopt_pop_heatmap_df_list <- list()

cat("\nPre-procesando datos crudos para todos los heatmaps (SD)...\n")
for (current_tau_sd_proc in TAU_NORMAL_SD_SWEEP_LIST) {
  sd_label_proc <- paste0("sd_", sprintf("%.2f", current_tau_sd_proc))
  current_sd_all_means_raw_results_proc <- all_sds_raw_results_list_from_file[[sd_label_proc]]
  if (is.null(current_sd_all_means_raw_results_proc)) next
  
  heatmap_data_tm_this_sd_list_proc <- list()
  heatmap_data_aa_this_sd_list_proc <- list()
  heatmap_data_ar_pop_this_sd_list_proc <- list()
  heatmap_data_as_pop_this_sd_list_proc <- list()
  
  for (current_threshold_mean_proc in THRESHOLD_MEAN_SWEEP_LIST) {
    mean_label_proc <- paste0("mean_", sprintf("%.2f", current_threshold_mean_proc))
    raw_df_this_mean_sd_proc <- current_sd_all_means_raw_results_proc[[mean_label_proc]]
    if (is.null(raw_df_this_mean_sd_proc) || nrow(raw_df_this_mean_sd_proc) == 0) next
    
    NUM_RUNS_THIS_COMBO_ACTUAL_PROC <- length(unique(raw_df_this_mean_sd_proc$run_id))
    if (NUM_RUNS_THIS_COMBO_ACTUAL_PROC == 0) next
    
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
      ungroup()
    
    panel_data_tm_list_proc <- list()
    panel_data_aa_list_proc <- list()
    panel_data_ar_pop_list_proc <- list()
    panel_data_as_pop_list_proc <- list()
    
    for (iul_val_proc in IUL_VALUES_SWEEP) {
      for (h_val_proc in H_VALUES_SWEEP) {
        current_cell_data_proc <- base_run_summary_proc %>%
          filter(innovation_iul_Gamma == iul_val_proc, social_distance_h == h_val_proc)
        
        if (nrow(current_cell_data_proc) == 0) next
        
        # --- SD 1: Adopción total
        sd_aa_cell <- sd(current_cell_data_proc$adopters_prop_at_cell, na.rm = TRUE)
        panel_data_aa_list_proc[[length(panel_data_aa_list_proc) + 1]] <- 
          data.frame(iul = iul_val_proc, h = h_val_proc, val = sd_aa_cell)
        
        # --- SD 2: Transición (variabilidad del indicador binario)
        sd_tm_cell <- sd(current_cell_data_proc$is_transition_vs_prev_gamma, na.rm = TRUE)
        panel_data_tm_list_proc[[length(panel_data_tm_list_proc) + 1]] <- 
          data.frame(iul = iul_val_proc, h = h_val_proc, val = sd_tm_cell)
        
        # --- SD 3: Adopción racional
        sd_ar_cell <- sd(current_cell_data_proc$prop_rational_of_pop_at_cell, na.rm = TRUE)
        panel_data_ar_pop_list_proc[[length(panel_data_ar_pop_list_proc) + 1]] <- 
          data.frame(iul = iul_val_proc, h = h_val_proc, val = sd_ar_cell)
        
        # --- SD 4: Adopción social
        sd_as_cell <- sd(current_cell_data_proc$prop_social_of_pop_at_cell, na.rm = TRUE)
        panel_data_as_pop_list_proc[[length(panel_data_as_pop_list_proc) + 1]] <- 
          data.frame(iul = iul_val_proc, h = h_val_proc, val = sd_as_cell)

        # Actualizar el máximo global de SD
        max_sd_global <- max(max_sd_global, sd_aa_cell, sd_tm_cell, sd_ar_cell, sd_as_cell, na.rm = TRUE)
      }
    }
    
    heatmap_data_tm_this_sd_list_proc[[mean_label_proc]] <-
      bind_rows(panel_data_tm_list_proc) %>%
      mutate(tau_mean_param = current_threshold_mean_proc,
             sd_transition_value_to_plot = val,
             tau_sd_param = current_tau_sd_proc) %>% select(-val)
    
    heatmap_data_aa_this_sd_list_proc[[mean_label_proc]] <-
      bind_rows(panel_data_aa_list_proc) %>%
      mutate(tau_mean_param = current_threshold_mean_proc,
             sd_adopters_prop_to_plot = val,
             tau_sd_param = current_tau_sd_proc) %>% select(-val)
    
    heatmap_data_ar_pop_this_sd_list_proc[[mean_label_proc]] <-
      bind_rows(panel_data_ar_pop_list_proc) %>%
      mutate(tau_mean_param = current_threshold_mean_proc,
             sd_rational_adopt_pop_to_plot = val,
             tau_sd_param = current_tau_sd_proc) %>% select(-val)
    
    heatmap_data_as_pop_this_sd_list_proc[[mean_label_proc]] <-
      bind_rows(panel_data_as_pop_list_proc) %>%
      mutate(tau_mean_param = current_threshold_mean_proc,
             sd_social_adopt_pop_to_plot = val,
             tau_sd_param = current_tau_sd_proc) %>% select(-val)
  }
  
  all_sds_transition_metric_heatmap_df_list[[sd_label_proc]] <-
    bind_rows(heatmap_data_tm_this_sd_list_proc)
  all_sds_avg_adoption_heatmap_df_list[[sd_label_proc]] <-
    bind_rows(heatmap_data_aa_this_sd_list_proc)
  all_sds_avg_rational_adopt_pop_heatmap_df_list[[sd_label_proc]] <-
    bind_rows(heatmap_data_ar_pop_this_sd_list_proc)
  all_sds_avg_social_adopt_pop_heatmap_df_list[[sd_label_proc]] <-
    bind_rows(heatmap_data_as_pop_this_sd_list_proc)
}
cat("Máximo global de SD detectado:", max_sd_global, "\n")
cat("Pre-procesamiento completado (SD por celda).\n")

# --- Función de creación de heatmap (idéntica a la versión v3 del original) ---
create_single_heatmap_v3 <- function(df_plot_data, fill_col_name, legend_title_text, viridis_option,
                                     show_legend=TRUE, y_axis_label_on=TRUE, x_axis_label_on=TRUE,
                                     panel_row_title="") {
  if (is.null(df_plot_data) || nrow(df_plot_data) == 0 || all(is.na(df_plot_data[[fill_col_name]]))) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No plottable data") +
             theme_void() +
             labs(y = if(y_axis_label_on) panel_row_title else NULL) +
             theme(axis.title.y = element_text(size=8, face="bold", angle=90)))
  }
  
  if ("iul" %in% names(df_plot_data)) df_plot_data <- rename(df_plot_data, innovation_iul_Gamma = iul)
  if ("h" %in% names(df_plot_data)) df_plot_data <- rename(df_plot_data, social_distance_h = h)
  
  h_levels_sorted <- sprintf("%.2f", sort(unique(H_VALUES_SWEEP)))
  df_plot_data$social_distance_h_factor <- factor(sprintf("%.2f", df_plot_data$social_distance_h), levels = h_levels_sorted)
  
  y_breaks <- sprintf("%.2f", H_VALUES_SWEEP[seq(1, length(H_VALUES_SWEEP), by=2)])
  y_labels <- y_breaks
  x_breaks <- seq(0, 1, 0.25)
  x_labels <- sprintf("%.2f", x_breaks)
  
  ggplot(df_plot_data, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = .data[[fill_col_name]])) +
    geom_tile(color="white", lwd=0.1) +
    scale_fill_viridis_c(
      name = if(show_legend) legend_title_text else NULL,
      limits = c(0, max_sd_global),  # Usar el límite global
      option = viridis_option,
      na.value = "grey90"
    ) +
    labs(x = if(x_axis_label_on) expression(paste("IUL (", Gamma, ")")) else NULL,
         y = if(y_axis_label_on) panel_row_title else NULL) +
    scale_x_continuous(breaks=x_breaks, labels=if(x_axis_label_on) x_labels else NULL, expand=c(0,0)) +
    scale_y_discrete(drop=FALSE, breaks=y_breaks, labels=if(y_axis_label_on) y_labels else NULL) +
    theme_minimal(base_size=7) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=6),
          axis.text.y=element_text(size=6),
          axis.title.x=element_text(size=8, face="bold"),
          axis.title.y=element_text(size=8, face="bold", angle=90),
          legend.position=if(show_legend) "right" else "none",
          legend.title=element_text(size=7),
          legend.text=element_text(size=6))
}

# --- Bucle principal de ploteo ---
for (current_threshold_mean_plot in THRESHOLD_MEAN_SWEEP_LIST) {
  cat(paste0("\nGenerating SD Heatmaps for Mean τ = ", current_threshold_mean_plot, "\n"))
  
  metric_titles_for_cols <- c("SD Adoption",
                              "SD Phase Trans.",
                              "SD Rational",
                              "SD Social")
  
  metric_fill_vars_for_cols <- c("sd_adopters_prop_to_plot",
                                 "sd_transition_value_to_plot",
                                 "sd_rational_adopt_pop_to_plot",
                                 "sd_social_adopt_pop_to_plot")
  
  metric_color_options_for_cols <- rep("magma", 4)
  
  data_sources_for_cols <- list(
    all_sds_avg_adoption_heatmap_df_list,
    all_sds_transition_metric_heatmap_df_list,
    all_sds_avg_rational_adopt_pop_heatmap_df_list,
    all_sds_avg_social_adopt_pop_heatmap_df_list
  )
  
  plot_list_for_this_mean_pdf_ordered <- list()
  
  for (row_idx in 1:length(TAU_NORMAL_SD_SWEEP_LIST)) {
    current_tau_sd_plot <- TAU_NORMAL_SD_SWEEP_LIST[row_idx]
    sd_label_plot <- paste0("sd_", sprintf("%.2f", current_tau_sd_plot))
    
    for (col_idx in 1:4) {
      df_all_means_for_sd_metric <- data_sources_for_cols[[col_idx]][[sd_label_plot]]
      current_df_for_panel <- if(!is.null(df_all_means_for_sd_metric)) {
        df_all_means_for_sd_metric %>% filter(tau_mean_param == current_threshold_mean_plot)
      } else data.frame(innovation_iul_Gamma=numeric(0), social_distance_h=numeric(0))
      
      row_title_str <- if(col_idx == 1) paste0("MSP (h) - SD=", sprintf("%.2f", current_tau_sd_plot)) else ""
      y_label_visible <- (col_idx == 1)
      x_label_visible <- (row_idx == length(TAU_NORMAL_SD_SWEEP_LIST))
      legend_visible <- (col_idx == 4)
      
      plot_index_in_list <- ((row_idx-1)*4) + col_idx
      
      plot_list_for_this_mean_pdf_ordered[[plot_index_in_list]] <-
        create_single_heatmap_v3(
          df_plot_data = current_df_for_panel,
          fill_col_name = metric_fill_vars_for_cols[col_idx],
          legend_title_text = NULL,
          viridis_option = metric_color_options_for_cols[col_idx],
          show_legend = legend_visible,
          y_axis_label_on = y_label_visible,
          x_axis_label_on = x_label_visible,
          panel_row_title = row_title_str
        )
    }
  }
  
  col_titles_plots <- lapply(metric_titles_for_cols, function(title) {
    ggplot() + labs(title=title) + theme_void() +
      theme(plot.title = element_text(hjust=0.5, size=10, face="bold"))
  })
  
  column_titles_row_layout <- Reduce(`+`, col_titles_plots) + plot_layout(ncol=4)
  heatmaps_grid_layout <- wrap_plots(plot_list_for_this_mean_pdf_ordered, ncol=4, byrow=TRUE)
  
  final_combined_layout <- column_titles_row_layout / heatmaps_grid_layout + plot_layout(heights=c(0.05, 1))
  
  final_plot_with_annotation <- final_combined_layout +
    plot_annotation(
      title = paste("SD Heatmaps for ATP-net - Mean threshold =", sprintf("%.2f", current_threshold_mean_plot)),
      subtitle = paste("Standard deviation across 96 runs per (IUL,h). Seeding strategy:", SEEDING_STRATEGY_FIXED),
      theme = theme(plot.title=element_text(hjust=0.5, face="bold", size=12),
                    plot.subtitle=element_text(hjust=0.5, size=9.5))
    )
  
  pdf_width <- 7.5
  pdf_height <- 7.0
  plot_filename_consolidated_final <- paste0(PLOTS_DIR, "heatmaps_SD__seed_",
                                             SEEDING_STRATEGY_FIXED, "_mean_tau_",
                                             sprintf("%.2f", current_threshold_mean_plot), ".pdf")
  ggsave(plot_filename_consolidated_final, final_plot_with_annotation,
         width=pdf_width, height=pdf_height, limitsize=FALSE)
  cat(paste0("  Saved SD PDF: ", plot_filename_consolidated_final, "\n"))
}

cat("\nGeneración de todos los PDFs SD completada.\n")