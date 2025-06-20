# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script para Analizar Resultados CRUDOS y Generar PDF Consolidados con 16 Heatmaps
# Heatmaps 3 y 4 AHORA muestran:
#   "Proporción Media de la POBLACIÓN Adoptada por Elección Racional"
#   "Proporción Media de la POBLACIÓN Adoptada por Influencia Social"
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

SEEDING_STRATEGY_FIXED <- "closeness"#"closeness" #"marginal" #"eigen" #"central" #"random"
PHASE_TRANSITION_THRESHOLD_JUMP <- 1/3 
SUCCESSFUL_DIFFUSION_THRESHOLD_PROP <- 0.50 
# MIN_PROP_SUCCESSFUL_RUNS_FOR_TILE_CELL ya no se usa para el heatmap de transiciones (heatmap 2)

cat("Cargando resultados crudos guardados...\n")
grand_raw_results_path <- paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_raw_results_all_sds_means_", SEEDING_STRATEGY_FIXED, ".rds") 
if (!file.exists(grand_raw_results_path)) {
  stop("Archivo de resultados crudos no encontrado: ", grand_raw_results_path)
}
all_sds_raw_results_list_from_file <- readRDS(grand_raw_results_path)
cat("Resultados cargados.\n")

# --- PASO 0: Pre-procesar los datos crudos ---
all_sds_transition_metric_heatmap_df_list <- list()
all_sds_avg_adoption_heatmap_df_list    <- list()
all_sds_avg_rational_adopt_pop_heatmap_df_list <- list() # NUEVO NOMBRE
all_sds_avg_social_adopt_pop_heatmap_df_list   <- list() # NUEVO NOMBRE

cat("\nPre-procesando datos crudos para todos los heatmaps...\n")
for (current_tau_sd_proc in TAU_NORMAL_SD_SWEEP_LIST) {
  sd_label_proc <- paste0("sd_", sprintf("%.2f", current_tau_sd_proc))
  current_sd_all_means_raw_results_proc <- all_sds_raw_results_list_from_file[[sd_label_proc]]
  if (is.null(current_sd_all_means_raw_results_proc)) { next }
  
  heatmap_data_tm_this_sd_list_proc <- list(); heatmap_data_aa_this_sd_list_proc <- list()
  heatmap_data_ar_pop_this_sd_list_proc <- list(); heatmap_data_as_pop_this_sd_list_proc <- list() # NUEVO
  
  for (current_threshold_mean_proc in THRESHOLD_MEAN_SWEEP_LIST) {
    mean_label_proc <- paste0("mean_", sprintf("%.2f", current_threshold_mean_proc))
    raw_df_this_mean_sd_proc <- current_sd_all_means_raw_results_proc[[mean_label_proc]]
    if (is.null(raw_df_this_mean_sd_proc) || nrow(raw_df_this_mean_sd_proc) == 0) { next }
    
    NUM_RUNS_THIS_COMBO_ACTUAL_PROC <- length(unique(raw_df_this_mean_sd_proc$run_id))
    if(NUM_RUNS_THIS_COMBO_ACTUAL_PROC == 0) next
    
    # Pre-cálculo general: adopters_prop, y ahora prop_rational_OF_POPULATION, prop_social_OF_POPULATION
    base_run_summary_proc <- raw_df_this_mean_sd_proc %>%
      group_by(run_id, social_distance_h, innovation_iul_Gamma) %>% 
      summarise(
        adopters_prop_at_cell = first(num_adopters) / first(N_nodes_actual),
        # NUEVAS MÉTRICAS: Proporción de la población total adoptada por cada mecanismo
        prop_rational_of_pop_at_cell = first(num_adopted_rational) / first(N_nodes_actual),
        prop_social_of_pop_at_cell = first(num_adopted_social) / first(N_nodes_actual),
        # Columnas necesarias para el heatmap de transición
        num_adopters_at_cell = first(num_adopters), # Necesaria para is_successful_at_cell si se usa
        n_nodes_at_cell = first(N_nodes_actual),    # Necesaria para is_successful_at_cell
        .groups = 'drop'
      ) %>%
      mutate(
        # Esta columna solo se usa para el filtro del heatmap de transiciones (si se reintroduce)
        # is_successful_at_cell = adopters_prop_at_cell >= SUCCESSFUL_DIFFUSION_THRESHOLD_PROP 
      ) %>%
      arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
      group_by(run_id, social_distance_h) %>%
      mutate(
        jump_at_step = adopters_prop_at_cell - lag(adopters_prop_at_cell), 
        is_transition_vs_prev_gamma = ifelse(!is.na(jump_at_step) & jump_at_step >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
      ) %>%
      group_by(run_id, social_distance_h) %>% # Re-agrupar para el siguiente mutate
      mutate(first_transition_IUL_for_series = if (any(is_transition_vs_prev_gamma == 1, na.rm=TRUE)) {min(innovation_iul_Gamma[which(is_transition_vs_prev_gamma == 1)])} else {NA_real_}) %>%
      ungroup() %>%
      # Ahora para las transiciones por H
      arrange(run_id, innovation_iul_Gamma, social_distance_h) %>%
      group_by(run_id, innovation_iul_Gamma) %>%
      mutate(
        jump_vs_prev_h = adopters_prop_at_cell - lag(adopters_prop_at_cell),
        is_transition_vs_prev_h = ifelse(!is.na(jump_vs_prev_h) & jump_vs_prev_h >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
      ) %>%
      group_by(run_id, innovation_iul_Gamma) %>% # Re-agrupar
      mutate(first_transition_H_for_series = if (any(is_transition_vs_prev_h == 1, na.rm=TRUE)) {min(social_distance_h[which(is_transition_vs_prev_h == 1)])} else {NA_real_}) %>%
      ungroup()
    
    
    panel_data_tm_list_proc <- list(); panel_data_aa_list_proc <- list()
    panel_data_ar_pop_list_proc <- list(); panel_data_as_pop_list_proc <- list() # NUEVO
    
    for (iul_val_proc in IUL_VALUES_SWEEP) {
      for (h_val_proc in H_VALUES_SWEEP) {
        
        # Metric 1: Transition Metric (NO ACUMULATIVA, sensible a Γ y h, SIN filtro de éxito de celda)
        runs_transitioned_to_this_cell_by_gamma <- base_run_summary_proc %>%
          filter(social_distance_h == h_val_proc, 
                 !is.na(first_transition_IUL_for_series), 
                 first_transition_IUL_for_series == iul_val_proc) %>% pull(run_id)
        runs_transitioned_to_this_cell_by_h <- base_run_summary_proc %>%
          filter(innovation_iul_Gamma == iul_val_proc, 
                 !is.na(first_transition_H_for_series), 
                 first_transition_H_for_series == h_val_proc) %>% pull(run_id)
        unique_runs_transitioned_to_cell <- unique(c(runs_transitioned_to_this_cell_by_gamma, runs_transitioned_to_this_cell_by_h))
        prop_tm_cell_combined <- length(unique_runs_transitioned_to_cell) / NUM_RUNS_THIS_COMBO_ACTUAL_PROC
        panel_data_tm_list_proc[[length(panel_data_tm_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_tm_cell_combined)
        
        # Métricas 2, 3 (NUEVO), 4 (NUEVO)
        current_cell_data_proc <- base_run_summary_proc %>%
          filter(innovation_iul_Gamma == iul_val_proc, social_distance_h == h_val_proc)
        
        # Metric 2: Average Adoption (total)
        prop_aa_cell <- NA_real_; if(nrow(current_cell_data_proc) == NUM_RUNS_THIS_COMBO_ACTUAL_PROC) {prop_aa_cell <- mean(current_cell_data_proc$adopters_prop_at_cell, na.rm=TRUE); if(is.nan(prop_aa_cell)) prop_aa_cell <- NA_real_}
        panel_data_aa_list_proc[[length(panel_data_aa_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_aa_cell)
        
        # Metric 3: Avg. Adoption by Rational Choice (prop. of population)
        prop_ar_pop_cell <- NA_real_; if(nrow(current_cell_data_proc) == NUM_RUNS_THIS_COMBO_ACTUAL_PROC) {prop_ar_pop_cell <- mean(current_cell_data_proc$prop_rational_of_pop_at_cell, na.rm=TRUE); if(is.nan(prop_ar_pop_cell)) prop_ar_pop_cell <- NA_real_}
        panel_data_ar_pop_list_proc[[length(panel_data_ar_pop_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_ar_pop_cell)
        
        # Metric 4: Avg. Adoption by Social Influence (prop. of population)
        prop_as_pop_cell <- NA_real_; if(nrow(current_cell_data_proc) == NUM_RUNS_THIS_COMBO_ACTUAL_PROC) {prop_as_pop_cell <- mean(current_cell_data_proc$prop_social_of_pop_at_cell, na.rm=TRUE); if(is.nan(prop_as_pop_cell)) prop_as_pop_cell <- NA_real_}
        panel_data_as_pop_list_proc[[length(panel_data_as_pop_list_proc) + 1]] <- data.frame(iul=iul_val_proc, h=h_val_proc, val=prop_as_pop_cell)
      }
    }
    heatmap_data_tm_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_tm_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, proportion_value_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)
    heatmap_data_aa_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_aa_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, mean_adopters_prop_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)
    heatmap_data_ar_pop_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_ar_pop_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, avg_rational_adopt_pop_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val) # NUEVO
    heatmap_data_as_pop_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_as_pop_list_proc) %>% mutate(tau_mean_param = current_threshold_mean_proc, avg_social_adopt_pop_to_plot = val, tau_sd_param = current_tau_sd_proc) %>% select(-val)   # NUEVO
  }
  all_sds_transition_metric_heatmap_df_list[[sd_label_proc]] <- bind_rows(heatmap_data_tm_this_sd_list_proc[!sapply(heatmap_data_tm_this_sd_list_proc, is.null)])
  all_sds_avg_adoption_heatmap_df_list[[sd_label_proc]]    <- bind_rows(heatmap_data_aa_this_sd_list_proc[!sapply(heatmap_data_aa_this_sd_list_proc, is.null)])
  all_sds_avg_rational_adopt_pop_heatmap_df_list[[sd_label_proc]] <- bind_rows(heatmap_data_ar_pop_this_sd_list_proc[!sapply(heatmap_data_ar_pop_this_sd_list_proc, is.null)]) # NUEVO
  all_sds_avg_social_adopt_pop_heatmap_df_list[[sd_label_proc]]   <- bind_rows(heatmap_data_as_pop_this_sd_list_proc[!sapply(heatmap_data_as_pop_this_sd_list_proc, is.null)])   # NUEVO
}
cat("Pre-procesamiento de todos los datos para heatmaps completado (con nuevas métricas de tipo de adopción).\n")

# --- Función Auxiliar para crear UN heatmap individual ---
create_single_heatmap_v3 <- function(df_plot_data, fill_col_name, legend_title_text, viridis_option, 
                                     show_legend=TRUE, y_axis_label_on=TRUE, x_axis_label_on=TRUE, 
                                     panel_row_title="") { # Cambiado panel_super_title a panel_row_title
  
  if (is.null(df_plot_data) || nrow(df_plot_data) == 0 || all(is.na(df_plot_data[[fill_col_name]]))) {
    # Si hay un título de fila, mostrarlo incluso para plot vacío, para mantener alineación
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No plottable data") + 
             theme_void() + 
             labs(title = NULL, y = if(y_axis_label_on) panel_row_title else NULL) + # Título de fila como etiqueta Y
             theme(axis.title.y = element_text(size=8, face="bold", angle=90, margin = margin(r=5))) # Estilo para el título de fila
    )
  }
  
  if("iul" %in% names(df_plot_data)) df_plot_data <- rename(df_plot_data, innovation_iul_Gamma = iul)
  if("h" %in% names(df_plot_data)) df_plot_data <- rename(df_plot_data, social_distance_h = h)
  
  # Asegurar que los niveles del factor h sean todos los posibles para consistencia del eje Y
  h_levels_sorted <- sprintf("%.2f", sort(unique(H_VALUES_SWEEP)))
  df_plot_data$social_distance_h_factor <- factor(sprintf("%.2f", df_plot_data$social_distance_h), levels = h_levels_sorted)
  
  # Definir breaks y labels para el eje Y (MSP h)
  y_breaks <- sprintf("%.2f", H_VALUES_SWEEP[seq(1, length(H_VALUES_SWEEP), by = 2)]) # Cada dos valores
  y_labels <- y_breaks
  
  # Definir breaks y labels para el eje X (IUL Gamma)
  x_breaks <- seq(0, 1, 0.25)
  x_labels <- sprintf("%.2f", x_breaks)
  
  plot_obj <- ggplot(df_plot_data, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = .data[[fill_col_name]])) +
    geom_tile(color = "white", lwd = 0.1) + # Líneas blancas más sutiles
    scale_fill_viridis_c(
      name = if(show_legend) legend_title_text else NULL, # Quitar título de leyenda si show_legend es FALSE o es NULL
      limits = c(0, 1), option=viridis_option, n.breaks=4, na.value = "grey90" # Celdas NA en gris claro
    ) +
    labs(
      x = if(x_axis_label_on) expression(paste("IUL (", Gamma, ")")) else NULL, 
      # Título de fila (SD) ahora se manejará con patchwork o como etiqueta Y del primer plot de la fila
      y = if(y_axis_label_on) panel_row_title else NULL, # Usar panel_row_title como etiqueta del eje Y
      title = NULL # Los títulos de columna se pondrán encima con patchwork
    ) +
    scale_x_continuous(breaks = x_breaks, labels= if(x_axis_label_on) x_labels else NULL, expand = c(0,0)) +
    scale_y_discrete(drop = FALSE, breaks = y_breaks, labels = if(y_axis_label_on) y_labels else NULL) + 
    theme_minimal(base_size = 7) + # Tamaño base general más pequeño
    theme(
      # plot.title = element_text(hjust = 0.5, face = "bold", size=8), # Ya no se usa aquí
      axis.text.x = element_text(angle = 45, hjust = 1, size=6, 
                                 color = if(x_axis_label_on) "black" else "transparent"), # Ocultar si no es la última fila
      axis.text.y = element_text(size=6, 
                                 color = if(y_axis_label_on) "black" else "transparent"), # Ocultar si no es la primera columna
      #axis.title.x = element_text(size=7, margin = margin(t = 2, unit="mm")), 
      axis.title.x = element_text(size=8, face="bold", margin = margin(t = 2, unit="mm")),
      axis.title.y = element_text(size=8, face="bold", angle=90, margin = margin(r = 2, unit="mm")), # Para el título de fila "SD = X.XX"
      legend.position = if(show_legend) "right" else "none",
      legend.title = element_text(size = 7), 
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.6, "lines"), 
      panel.grid = element_line(linewidth = 0.15), 
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")
    )
  return(plot_obj)
}


# --- Bucle Principal de PLOTEO: Generar un PDF por cada MEAN_tau ---
for (current_threshold_mean_plot in THRESHOLD_MEAN_SWEEP_LIST) {
  cat(paste0("\nGenerating Consolidated PDF for Mean τ = ", current_threshold_mean_plot, "\n"))
  
  plot_list_for_this_mean_pdf_ordered <- list() 
  
  # Definiciones de las columnas (Mismo orden que en tu PDF de ejemplo)
  metric_titles_for_cols <- c("Avg. Adoption", 
                              "Phase Trans Prob.",
                              "Avg. Adopt.\nby Rational", 
                              "Avg. Adopt.\nby Social Infl.")
  
  metric_fill_vars_for_cols <- c("mean_adopters_prop_to_plot",    
                                 "proportion_value_to_plot",      
                                 "avg_rational_adopt_pop_to_plot",
                                 "avg_social_adopt_pop_to_plot")  
  
  metric_legend_titles_for_cols <- c("Mean Total\nAdopt. Prop.", # Leyenda para Col 1
                                     "Prop. Runs with\n1st Transition", # Leyenda para Col 2
                                     "Mean Pop. Adopted\nvia Rational",  # Leyenda para Col 3
                                     "Mean Pop. Adopted\nvia Social Infl.") # Leyenda para Col 4
  
  metric_color_options_for_cols <- rep("viridis", 4) # Punto 3: Todos viridis
  
  # Orden de las fuentes de datos debe coincidir con el orden de las columnas arriba
  data_sources_for_cols <- list(
    all_sds_avg_adoption_heatmap_df_list,       # Para "Avg. Adoption"
    all_sds_transition_metric_heatmap_df_list,  # Para "Phase Trans Prob."
    all_sds_avg_rational_adopt_pop_heatmap_df_list, # Para "Avg. Adopt. by Rational"
    all_sds_avg_social_adopt_pop_heatmap_df_list    # Para "Avg. Adopt. by Social Infl."
  )
  
  # Bucle sobre las Desviaciones Estándar (Filas del PDF)
  for (row_idx in 1:length(TAU_NORMAL_SD_SWEEP_LIST)) {
    current_tau_sd_plot <- TAU_NORMAL_SD_SWEEP_LIST[row_idx]
    sd_label_plot <- paste0("sd_", sprintf("%.2f", current_tau_sd_plot))
    # cat(paste0("  Preparing row for SD τ = ", current_tau_sd_plot, "...\n")) # Verboso
    
    # Crear los 4 heatmaps para esta fila (este SD_tau)
    for (col_idx in 1:4) { 
      # Extraer el dataframe de datos pre-procesados para el SD y la Métrica actual
      df_all_means_for_sd_metric <- data_sources_for_cols[[col_idx]][[sd_label_plot]]
      
      # Filtrar por la media actual del PDF
      current_df_for_panel <- if(!is.null(df_all_means_for_sd_metric) && nrow(df_all_means_for_sd_metric) > 0) {
        df_all_means_for_sd_metric %>% filter(tau_mean_param == current_threshold_mean_plot)
      } else { 
        # Crear un dataframe vacío con las columnas esperadas si no hay datos
        # para que create_single_heatmap_v3 pueda manejarlo y mostrar "No plottable data"
        data.frame(innovation_iul_Gamma=numeric(0), social_distance_h=numeric(0)) %>% 
          mutate(!!metric_fill_vars_for_cols[col_idx] := numeric(0) ) # Añadir la columna de fill
      }
      
      # Título de Fila (SD) solo para el primer plot de la fila (col_idx == 1)
      #row_title_str <- if(col_idx == 1) paste0("SD τ=",sprintf("%.2f", current_tau_sd_plot)) else ""
      row_title_str <- if(col_idx == 1) paste0("MSP (h) - SD=",sprintf("%.2f", current_tau_sd_plot)) else ""
      
      y_label_visible <- (col_idx == 1) # Etiqueta Y ("MSP (h)" + SD) solo para la primera columna
      x_label_visible <- (row_idx == length(TAU_NORMAL_SD_SWEEP_LIST)) # Etiqueta X solo para la última fila
      legend_visible <- (col_idx == 4) # Leyenda solo para la última columna de cada fila
      
      plot_index_in_list <- ( (row_idx-1)*4 ) + col_idx
      
      plot_list_for_this_mean_pdf_ordered[[plot_index_in_list]] <- 
        create_single_heatmap_v3( # Usando v3
          df_plot_data = current_df_for_panel, 
          fill_col_name = metric_fill_vars_for_cols[col_idx], 
          legend_title_text = NULL, # Punto 4: Quitar título de la leyenda individual del heatmap
          viridis_option = metric_color_options_for_cols[col_idx],
          show_legend = legend_visible,
          y_axis_label_on = y_label_visible,
          x_axis_label_on = x_label_visible,
          panel_row_title = row_title_str # Este será el título de la fila, usado como etiqueta Y
        )
    } # Fin bucle metricas (columnas)
  } # Fin bucle SD_tau (filas)
  
  # --- Ensamblaje con Patchwork ---
  if (length(plot_list_for_this_mean_pdf_ordered) == (length(TAU_NORMAL_SD_SWEEP_LIST) * 4) ) {
    
    col_titles_plots <- lapply(metric_titles_for_cols, function(title) {
      ggplot() + labs(title=title) + theme_void() + 
        theme(plot.title = element_text(hjust=0.5, size=10, face="bold", margin = margin(b=0, t=2, unit="mm"))) # Tamaño y margen título columna
    })
    
    column_titles_row_layout <- Reduce(`+`, col_titles_plots) + plot_layout(ncol = 4)
    heatmaps_grid_layout <- wrap_plots(plot_list_for_this_mean_pdf_ordered, ncol = 4, byrow = TRUE)
    
    # Ajustar alturas relativas: más para los plots, menos para los títulos de columna
    final_combined_layout <- column_titles_row_layout / heatmaps_grid_layout + 
      plot_layout(heights = c(0.05, 1)) # Fila de títulos más pequeña
    
    # Obtener num_runs para subtítulo
    num_runs_val_subtitle <- NA 
    first_valid_sd_label <- names(all_sds_raw_results_list_from_file)[which(!sapply(all_sds_raw_results_list_from_file, is.null))[1]]
    if(!is.na(first_valid_sd_label) && !is.null(all_sds_raw_results_list_from_file[[first_valid_sd_label]])){
      first_valid_mean_label_for_runs <- names(all_sds_raw_results_list_from_file[[first_valid_sd_label]])[which(!sapply(all_sds_raw_results_list_from_file[[first_valid_sd_label]], is.null))[1]]
      if(!is.na(first_valid_mean_label_for_runs) && !is.null(all_sds_raw_results_list_from_file[[first_valid_sd_label]][[first_valid_mean_label_for_runs]])){
        num_runs_val_subtitle <- length(unique(all_sds_raw_results_list_from_file[[first_valid_sd_label]][[first_valid_mean_label_for_runs]]$run_id))
      }
    }
    if(is.na(num_runs_val_subtitle)) num_runs_val_subtitle <- "N/A" # Fallback si no se encuentran datos
    
    # Añadir título y subtítulo general al PDF
    final_plot_with_annotation <- final_combined_layout + 
      plot_annotation(
        title = paste("Consolidated Heatmaps for ATP-net - Mean threshold =", sprintf("%.2f", current_threshold_mean_plot)),
        subtitle = paste("Thresholds ~ N(mu=", sprintf("%.2f", current_threshold_mean_plot), ", SD=var). ", 
                         num_runs_val_subtitle, " runs per (IUL,h) per individual panel. Seeding strategy: ", SEEDING_STRATEGY_FIXED,
                         sep=""),
        theme = theme(plot.title = element_text(hjust = 0.5, face="bold", size=12), # Tamaño título principal ajustado
                      plot.subtitle = element_text(hjust = 0.5, size=9.5)) # Tamaño subtítulo ajustado
      ) 
    
    pdf_width <- 7.5  # Mantenido de tu ejemplo
    pdf_height <- 7.0 # Ligeramente ajustado para el espaciado
    
    plot_filename_consolidated_final <- paste0(PLOTS_DIR, "heatmaps_mean_tau_", sprintf("%.2f", current_threshold_mean_plot), "_seed_",SEEDING_STRATEGY_FIXED, ".pdf")
    ggsave(plot_filename_consolidated_final, final_plot_with_annotation, width = pdf_width, height = pdf_height, limitsize = FALSE)
    cat(paste0("  Saved consolidated PDF: ", plot_filename_consolidated_final, "\n"))
    
  } else {
    cat(paste0("  No se generaron suficientes plots para el PDF consolidado de Mean τ = ", current_threshold_mean_plot, ".\n"))
  }
} 
cat("\nGeneración de todos los PDFs consolidados completada.\n")
 
# Guardar los datos pre-plot finales
# all_sds_transition_metric_heatmap_df_list no cambió su LÓGICA de cálculo, solo que ya no se filtra por éxito de celda
saveRDS(all_sds_transition_metric_heatmap_df_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_TM_NO_CELL_FILTER_heatmap_data_all_sds.rds"))
saveRDS(all_sds_avg_adoption_heatmap_df_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_AVG_ADOPTION_heatmap_data_all_sds_v3.rds")) # v3 para indicar que es de este run
saveRDS(all_sds_avg_rational_adopt_pop_heatmap_df_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_AVG_RATIONAL_POP_heatmap_data_all_sds.rds"))
saveRDS(all_sds_avg_social_adopt_pop_heatmap_df_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_AVG_SOCIAL_POP_heatmap_data_all_sds.rds"))

cat("Datos para heatmaps (New Rational/Social metrics) guardados.\n")