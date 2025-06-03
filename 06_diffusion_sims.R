# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Descripción: Script principal para ejecutar simulaciones de contagio complejo.
#              Llama a funciones definidas en 'simulation_functions.R'.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(igraph)
library(doParallel)
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork) # Para el plot final
library(intergraph)
library(cluster)

# Cargar el script con las funciones (asegúrate de que esté en el mismo directorio o ajusta la ruta)
source("diffusion_tests/simulation_functions.R") 

# -----------------------------------------------------------------------------
# 1. Generación e Importación de Topologías de Red
# -----------------------------------------------------------------------------

ATP_NETWORK <- TRUE

if (ATP_NETWORK) { # Carga de redes Small-World SDA desde archivos
  networks_dir <- "trabajo_1_files/ATP_network_ergm/"
  graphs_ATP <- list()

  for (i in 1:5) {
    ATP_net <- readRDS(paste0(networks_dir, "ATP_network_simulated_1000_mur_", sprintf("%03d", i), ".rds"))
    ATP_net <- asIgraph(ATP_net)
    
    graphs_ATP[[i]] <- ATP_net
  }
  
  graphs_list_to_simulate <- graphs_ATP
  current_graph_type_label <- "ATP-net"
  
} else { # Generación de redes sintéticas
  
  N_nodes_sim <- 1000
  density_sim <- 0.029
  num_edges_sim <- round(density_sim * N_nodes_sim * (N_nodes_sim - 1) / 2)
  m_sim <- 5          # Para Barabasi-Albert
  p_rewire_sim <- 0.2 # Para Small-World
  
  # graphs_list_to_simulate <- generate_small_world_networks(N_nodes_sim, num_edges_sim, p_rewire_sim)
  # current_graph_type_label <- "Small-World"
  
  graphs_list_to_simulate <- generate_scale_free_networks(N_nodes_sim, m_sim)
  current_graph_type_label <- "Scale-Free"
  
  # graphs_list_to_simulate <- generate_erdos_renyi_networks(N_nodes_sim, num_edges_sim)
  # current_graph_type_label <- "Erdos-Renyi"
}

# -----------------------------------------------------------------------------
# 2. Definición de Parámetros Globales para la Simulación
# -----------------------------------------------------------------------------
N_nodes_global <- vcount(graphs_list_to_simulate[[1]])

# Parámetros del espacio de simulación
homoph_values_sim <- seq(0.0, 0.3, length.out = 6) 
IUL_values_sim  <- seq(0.0, 1.0, by = 0.01) # Reducido para testeo más rápido
num_adopters_min_sim <- 0.1 # Proporción mínima para considerar éxito

# Thresholds
threshold_mean_list <- c(0.25, 0.30, 0.35) 
TAU_NORMAL_DISTRIBUTION_SD <- 0.16 

# Estrategia de selección de semillas: "PLci_top" o "random"
SEEDING_STRATEGY <- "random" # o "PLci_top"
NUM_INITIAL_SEEDS <- 8 # Cuántos de los top PLci probar como semilla
NUM_SUCCESSFUL_SIMS_PER_GRAPH <- 1 # Cuántas simulaciones "exitosas" guardar por grafo/umbral

# -----------------------------------------------------------------------------
# 3. Bucle Principal de Simulación
# -----------------------------------------------------------------------------
       
all_simulation_results_collection <- list() 

cat(paste("Iniciando simulaciones para tipo de red:", current_graph_type_label, "\n"))

for (current_threshold_mean in threshold_mean_list) { # τ base (fraccional 0-1)
  # T_str_sim <- as.character(current_threshold_mean)
  # cat(paste("  Corriendo para Umbral Base τ (fraccional):", T_str_sim, "\n"))
  T_str_sim <- paste0("norm_mean", sprintf("%.2f", current_threshold_mean), "_sd", sprintf("%.2f", TAU_NORMAL_DISTRIBUTION_SD))
  cat(paste("  Corriendo para Umbral τ ~ Normal(μ=", sprintf("%.2f", current_threshold_mean), ", σ=", TAU_NORMAL_DISTRIBUTION_SD, ")\n"))
  
  results_for_this_base_tau_and_all_graphs <- list() 
  
  for (graph_idx in seq_along(graphs_list_to_simulate)) {
    current_graph_obj_sim <- graphs_list_to_simulate[[graph_idx]]
    cat(paste("    Procesando grafo #", graph_idx, "de", length(graphs_list_to_simulate), "...\n"))
    
    node_mur_q_for_sim <- V(current_graph_obj_sim)$q_i # q_i (MUR) [CAMBIAR ALPHA --> Q_I EN LAS REDES SINTÉTICAS !!!!!]
    node_degrees_for_sim <- igraph::degree(current_graph_obj_sim)
    
    # Calcular la matriz de disimilitud de Gower
    attributes_for_distance <- data.frame(
      age = V(current_graph_obj_sim)$age,
      educ_num = V(current_graph_obj_sim)$educ_num,
      race = as.factor(V(current_graph_obj_sim)$race), # de 'character' a 'factor'
      relig = as.factor(V(current_graph_obj_sim)$relig),
      sex = as.factor(V(current_graph_obj_sim)$sex)
    )
    social_distance_matrix_d_ij <- as.matrix(daisy(attributes_for_distance, metric = "gower"))
    
    
    
    # Umbrales τ_i usando la distribución Normal
    set.seed(123 + graph_idx) 
    
    node_individual_thresholds_tau_frac_for_sim <- rnorm(
      n = N_nodes_global, # Número de nodos del grafo actual
      mean = current_threshold_mean,
      sd = TAU_NORMAL_DISTRIBUTION_SD
    )
    
    node_individual_thresholds_tau_frac_for_sim[node_individual_thresholds_tau_frac_for_sim < 0] <- 0
    node_individual_thresholds_tau_frac_for_sim[node_individual_thresholds_tau_frac_for_sim > 1] <- 1
    
    
    
    
    # Calcular umbrales de CONTEO para PLCI y para determinar el tamaño del clúster inicial
    # Estos son τ_conteo = τ_fraccional * grado_nodo (redondeado)
    node_thresholds_count_for_plci_and_cluster <- round(node_individual_thresholds_tau_frac_for_sim * node_degrees_for_sim)
    node_thresholds_count_for_plci_and_cluster[node_thresholds_count_for_plci_and_cluster <= 0] <- 1 # Mínimo 1
    
    # Selección de qué nodos semilla principales se van a probar
    seed_nodes_to_test_as_primary <- c()
    if (SEEDING_STRATEGY == "PLci_top") {
      adj_mat_for_plci_calc <- as.matrix(as_adjacency_matrix(current_graph_obj_sim))
      plci_values <- sapply(1:N_nodes_global, function(node_idx_for_plci) {
        get_PLci(
          seed_node = node_idx_for_plci, N_nodes = N_nodes_global, graph_obj = current_graph_obj_sim, 
          adj_matrix = adj_mat_for_plci_calc, 
          node_thresholds = node_thresholds_count_for_plci_and_cluster, # Umbral de CONTEO para PLCI
          num_initial_cluster_seeds_for_node = node_thresholds_count_for_plci_and_cluster[node_idx_for_plci], # Tamaño CONTEO para PLCI
          model_output_list_placeholder = NULL
        )[5] 
      })
      seed_nodes_to_test_as_primary <- order(as.numeric(plci_values), decreasing = TRUE)[1:min(N_nodes_global, NUM_INITIAL_SEEDS)]
      cat(paste("      Nodos semilla principales (PLci_top):", paste(seed_nodes_to_test_as_primary, collapse=", "), "\n"))
    } else if (SEEDING_STRATEGY == "random") {
      set.seed(graph_idx * sum(as.numeric(charToRaw(T_str_sim))) * 7) 
      seed_nodes_to_test_as_primary <- as.integer(sample(V(current_graph_obj_sim), NUM_INITIAL_SEEDS, replace = FALSE))
      cat(paste("      Nodos semilla principales (random):", paste(seed_nodes_to_test_as_primary, collapse=", "), "\n"))
    }
    
    if (length(seed_nodes_to_test_as_primary) == 0) {
      cat("      No se seleccionaron nodos semilla principales para este grafo. Saltando.\n")
      next
    }
    
    cl <- makeCluster(NUM_INITIAL_SEEDS, type = "FORK") # makeCluster(8, type = "FORK") M4
    registerDoParallel(cl)
    
    list_of_dfs_from_parallel_seeds <- foreach(
      current_primary_seed_id = seed_nodes_to_test_as_primary, 
      .combine = 'list',
      .multicombine = TRUE,
      .packages = c('igraph', 'dplyr'), 
      .export = c('sweep_homoph_parameter', 'get_complex_plot', 
                  'N_nodes_global', 'current_graph_obj_sim', 
                  'node_individual_thresholds_tau_frac_for_sim', # τ fraccional para get_complex_plot
                  'node_thresholds_count_for_plci_and_cluster', # τ conteo para clúster inicial
                  'node_mur_q_for_sim', 'node_degrees_for_sim',
                  'IUL_values_sim', 'homoph_values_sim', 'social_distance_matrix_d_ij'),
      .errorhandling = 'pass' 
    ) %dopar% {
      
      # Determinar el vector de semillas iniciales para ESTA simulación concreta
      # El tamaño del clúster inicial será el umbral de CONTEO del nodo semilla principal.
      # Esto asegura que la semilla principal "recibe" suficiente refuerzo inicial.
      num_seeds_for_initial_cluster = node_thresholds_count_for_plci_and_cluster[current_primary_seed_id]
      # Asegurar que no sea mayor que el número de nodos o grado + 1
      num_seeds_for_initial_cluster = min(num_seeds_for_initial_cluster, N_nodes_global, node_degrees_for_sim[current_primary_seed_id] + 1)
      if(num_seeds_for_initial_cluster < 1) num_seeds_for_initial_cluster <- 1

      
      initial_infectors_for_this_run <- c(current_primary_seed_id) # Empezar con la semilla principal
      
      if (num_seeds_for_initial_cluster > 1) {
        neighbors_of_primary <- as.numeric(neighbors(current_graph_obj_sim, current_primary_seed_id, mode="total"))
        # Necesitamos seleccionar (num_seeds_for_initial_cluster - 1) vecinos adicionales
        num_additional_needed = num_seeds_for_initial_cluster - 1
        
        if (length(neighbors_of_primary) >= num_additional_needed) {
          initial_infectors_for_this_run <- c(initial_infectors_for_this_run, sample(neighbors_of_primary, num_additional_needed))
        } else { # No hay suficientes vecinos, tomar todos y luego aleatorios si aún faltan
          initial_infectors_for_this_run <- c(initial_infectors_for_this_run, neighbors_of_primary)
          still_needed_more <- num_seeds_for_initial_cluster - length(initial_infectors_for_this_run)
          if (still_needed_more > 0) {
            potential_others <- setdiff(1:N_nodes_global, initial_infectors_for_this_run)
            if (length(potential_others) > 0) {
              initial_infectors_for_this_run <- c(initial_infectors_for_this_run, sample(potential_others, min(length(potential_others), still_needed_more)))
            }
          }
        }
      }
      initial_infectors_for_this_run <- unique(initial_infectors_for_this_run)
      
      df_from_worker <- sweep_homoph_parameter(
        primary_seed_id_arg = current_primary_seed_id,
        N_nodes_arg = N_nodes_global,
        graph_obj_arg = current_graph_obj_sim,
        node_individual_thresholds_tau_arg = node_individual_thresholds_tau_frac_for_sim, # τ fraccional para get_complex_plot
        node_mur_q_arg = node_mur_q_for_sim,
        all_innovation_iul_Gamma_values = IUL_values_sim, 
        all_social_distance_h_values = homoph_values_sim,   
        initial_infectors_vector_arg = initial_infectors_for_this_run, # El clúster inicial de nodos
        d_ij_matrix = social_distance_matrix_d_ij
      )
      
      df_from_worker$graph_type <- current_graph_type_label
      df_from_worker$base_threshold_tau_frac <- current_threshold_mean 
      df_from_worker$graph_instance_idx <- graph_idx
      
      return(df_from_worker)
    } 
    
    stopCluster(cl)
    
    final_dfs_to_combine_for_this_graph_instance <- list()
    successful_runs_count_for_filter = 0
    
    for(df_from_a_seed_run_parallel in list_of_dfs_from_parallel_seeds){
      if (inherits(df_from_a_seed_run_parallel, "simpleError")) {
        next 
      }
      #if (is.null(df_from_a_seed_run_parallel) || nrow(df_from_a_seed_run_parallel) == 0) next
      if (is.null(df_from_a_seed_run_parallel)) next
      
      process_this_df = FALSE
      if (SEEDING_STRATEGY == "PLci_top") {
        if (successful_runs_count_for_filter < NUM_SUCCESSFUL_SIMS_PER_GRAPH) {
          if (any(df_from_a_seed_run_parallel$num_adopters / N_nodes_global > num_adopters_min_sim, na.rm=TRUE)) {
            process_this_df = TRUE
          }
        }
      } else { 
        if (successful_runs_count_for_filter < NUM_SUCCESSFUL_SIMS_PER_GRAPH) {
          process_this_df = TRUE
        }
      }
      
      if(process_this_df){
        final_dfs_to_combine_for_this_graph_instance[[length(final_dfs_to_combine_for_this_graph_instance) + 1]] <- df_from_a_seed_run_parallel
        successful_runs_count_for_filter <- successful_runs_count_for_filter + 1
      }
    }
    
    if (length(final_dfs_to_combine_for_this_graph_instance) > 0) {
      results_for_this_base_tau_and_all_graphs[[paste0("graph_", graph_idx)]] <- bind_rows(final_dfs_to_combine_for_this_graph_instance)
    } else {
      cat(paste("    No se obtuvieron resultados válidos para el grafo #", graph_idx, "y umbral τ", T_str_sim, "\n"))
    }
  } 
  
  if (length(results_for_this_base_tau_and_all_graphs) > 0) {
    all_simulation_results_collection[[T_str_sim]] <- bind_rows(results_for_this_base_tau_and_all_graphs)
  }
} 

# Guardamos
saveRDS(all_simulation_results_collection, "trabajo_1_files/diffusion_simulation_files/diffusion_ATP_nets_1.rds")
all_simulation_results_collection <- readRDS("trabajo_1_files/diffusion_simulation_files/diffusion_ATP_nets_Random_h03-06_BIG.rds")

# -----------------------------------------------------------------------------
# 4. Visualización Rápida de Resultados 
# -----------------------------------------------------------------------------
if (length(all_simulation_results_collection) > 0) {
  plots_per_tau_threshold_comparison <- list()
  
  for (tau_label_for_plot_str in names(all_simulation_results_collection)) {
    data_for_this_tau_plot <- all_simulation_results_collection[[tau_label_for_plot_str]]
    
    if (!is.null(data_for_this_tau_plot) && nrow(data_for_this_tau_plot) > 0 &&
        all(c("innovation_iul_Gamma", "num_adopters", "seed", "social_distance_h", "graph_instance_idx") %in% names(data_for_this_tau_plot))) {
      
      plot_title_text <- paste(current_graph_type_label, "- Umbral Base τ (frac):", tau_label_for_plot_str)
      
      # Asegurar que social_distance_h sea factor para colores discretos
      # Y que los niveles estén ordenados para que la leyenda y colores coincidan
      h_levels <- unique(sort(data_for_this_tau_plot$social_distance_h))
      data_for_this_tau_plot$social_distance_h_factor <- factor(
        sprintf("%.2f", data_for_this_tau_plot$social_distance_h),
        levels = sprintf("%.2f", h_levels) 
      )
      h_labels_for_legend <- sprintf("%.2f", h_levels)
      
      fig_single_tau_comparison <- ggplot(
        data = data_for_this_tau_plot, 
        aes(x = innovation_iul_Gamma, 
            y = num_adopters / N_nodes_global, # Proporción de adoptadores
            color = social_distance_h_factor, 
            group = interaction(seed, social_distance_h_factor, graph_instance_idx) # Agrupar por corrida
        )
      ) +
        geom_line(linewidth = 0.4, alpha = 0.5) + # Líneas para cada corrida individual
        geom_point(size = 1.0, alpha = 0.5) +    # Puntos para cada dato individual
        
        # Paleta de colores Viridis "plasma", con dirección invertida para que amarillo sea h=0.00
        # El valor original de 'direction = -1' hacía que 0.00 fuera el color más oscuro (morado/azul).
        # 'direction = 1' debería invertirlo.
        # Ajustar 'begin' y 'end' puede ser necesario para afinar los colores exactos.
        scale_color_viridis_d(
          option = "plasma", 
          name = "Distancia Social (h)", 
          labels = h_labels_for_legend,
          direction = -1,  # Invertir la dirección de la paleta
          begin = 0.0,    # Empezar desde el inicio de la paleta (amarillo en plasma)
          end = 0.9       # Terminar un poco antes para evitar el morado más oscuro si se desea
        ) +
        labs(
          x = expression(paste("Utilidad Intrínseca Innovación (", Gamma, ")")), 
          y = "Proporción de Adoptadores", # Ya no es "Media"
          title = plot_title_text
        ) +
        ylim(0,1) + 
        theme_minimal(base_size = 10) +
        theme(
          legend.position = "bottom", 
          legend.direction = "horizontal", 
          plot.title = element_text(hjust = 0.5, size=11)
        )
      
      plots_per_tau_threshold_comparison[[plot_title_text]] <- fig_single_tau_comparison
      print(fig_single_tau_comparison) 
      
    } else {
      cat(paste("No hay datos suficientes o faltan columnas para graficar el umbral τ:", tau_label_for_plot_str, "\n"))
      if(!is.null(data_for_this_tau_plot)) print(paste("Nombres de columnas:", paste(names(data_for_this_tau_plot), collapse=", ")))
    }
  }
  
  # Combinar los gráficos para el tipo de red actual (si hay más de un umbral τ)
  if (length(plots_per_tau_threshold_comparison) > 0) {
    final_combined_plot_for_comparison <- wrap_plots(plots_per_tau_threshold_comparison, ncol = 1, guides = "collect") &
      theme(legend.position = "bottom")
    print(final_combined_plot_for_comparison)
    # Podrías guardar este gráfico combinado también:
    # filename_comparison_plot <- paste0("comparison_plot_", gsub("[:/ ]", "_", Sys.time()), "_", current_graph_type_label, ".pdf")
    # ggsave(filename_comparison_plot, final_combined_plot_for_comparison, width=7, height=2.5*length(plots_per_tau_threshold_comparison), units="in")
    cat(paste("Gráfico de comparación generado para:", current_graph_type_label, "\n"))
  }
  
} else {
  cat("No se generaron resultados consolidados para visualizar.\n")
}
