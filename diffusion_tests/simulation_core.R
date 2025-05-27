# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script 1: simulation_core.R
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

# Cargar el script con las funciones (asegúrate de que esté en el mismo directorio o ajusta la ruta)
source("diffusion_tests/simulation_functions.R") 

# -----------------------------------------------------------------------------
# 1. Generación e Importación de Topologías de Red
# -----------------------------------------------------------------------------

USE_SDA_FROM_FILES <- FALSE

if (USE_SDA_FROM_FILES) { # Carga de redes Small-World SDA desde archivos
  networks_dir <- "Talaga-homophily-network/"
  graphs_sda <- list()
  attributes_sda <- list()
  
  for (i in 1:5) {
    edge_file <- paste0(networks_dir, "talaga-homophily-network-edges-N324_", i, ".csv")
    edges <- read.csv(edge_file)
    attribute_file <- paste0(networks_dir, "talaga-homophily-network-node-attributes-N324_", i, ".csv")
    node_attributes <- read.csv(attribute_file)
    
    g <- graph_from_data_frame(d = edges, directed = FALSE, vertices = node_attributes)
    V(g)$alpha <- node_attributes$relevant_dim # REVISAR si 'relevant_dim' es 'alpha'
    
    graphs_sda[[i]] <- g
    attributes_sda[[i]] <- node_attributes # Guardar atributos
  }
  
  graphs_list_to_simulate <- graphs_sda
  current_graph_type_label <- "Small-World-SDA"
  
} else { # Generación de redes sintéticas
  
  N_nodes_sim <- 324
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
alpha_values_sim  <- seq(0.0, 1.0, by = 0.04) # Reducido para testeo más rápido
num_adopters_min_sim <- 0.1 # Proporción mínima para considerar éxito

# Tipo de umbral y distribución
T_dist_sim <- "homo"  # "homo" o "hetero"
T_type_sim <- "frac"  # "abs" o "frac"

threshold_values_list_sim <- c(0.10, 0.15, 0.20) 

# Estrategia de selección de semillas: "PLci_top" o "random"
SEEDING_STRATEGY <- "PLci_top" # o "random"
NUM_INITIAL_SEEDS_RANDOM <- 5 # Si SEEDING_STRATEGY es "random"
NUM_TOP_PLCI_NODES_TO_TEST <- 5 # Cuántos de los top PLci probar como semilla
NUM_SUCCESSFUL_SIMS_PER_GRAPH <- 1 # Cuántas simulaciones "exitosas" guardar por grafo/umbral

# -----------------------------------------------------------------------------
# 3. Bucle Principal de Simulación
# -----------------------------------------------------------------------------

all_simulation_results_collection <- list() 

cat(paste("Iniciando simulaciones para tipo de red:", current_graph_type_label, "\n"))

for (current_threshold_base_tau_fractional in threshold_values_list_sim) { # τ base (fraccional 0-1)
  T_str_sim <- as.character(current_threshold_base_tau_fractional)
  cat(paste("  Corriendo para Umbral Base τ (fraccional):", T_str_sim, "\n"))
  
  results_for_this_base_tau_and_all_graphs <- list() 
  
  for (graph_idx in seq_along(graphs_list_to_simulate)) {
    current_graph_obj_sim <- graphs_list_to_simulate[[graph_idx]]
    cat(paste("    Procesando grafo #", graph_idx, "de", length(graphs_list_to_simulate), "...\n"))
    
    node_mur_q_for_sim <- V(current_graph_obj_sim)$alpha # q_i (MUR)
    node_degrees_for_sim <- degree(current_graph_obj_sim, mode = "total")
    
    # Calcular umbrales individuales τ_i
    node_individual_thresholds_tau_frac_for_sim <- rep(current_threshold_base_tau_fractional, N_nodes_global)
    # (Aquí iría lógica para τ_i heterogéneos si T_dist_sim == "hetero")
    node_individual_thresholds_tau_frac_for_sim[node_individual_thresholds_tau_frac_for_sim <= 0 & current_threshold_base_tau_fractional > 0] <- 1e-6
    node_individual_thresholds_tau_frac_for_sim[node_individual_thresholds_tau_frac_for_sim > 1] <- 1.0
    
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
      seed_nodes_to_test_as_primary <- order(as.numeric(plci_values), decreasing = TRUE)[1:min(N_nodes_global, NUM_TOP_PLCI_NODES_TO_TEST)]
      cat(paste("      Nodos semilla principales (PLci_top):", paste(seed_nodes_to_test_as_primary, collapse=", "), "\n"))
    } else if (SEEDING_STRATEGY == "random") {
      set.seed(graph_idx * sum(as.numeric(charToRaw(T_str_sim))) * 7) 
      num_random_runs_per_graph = NUM_SUCCESSFUL_SIMS_PER_GRAPH 
      seed_nodes_to_test_as_primary <- sample(V(current_graph_obj_sim), num_random_runs_per_graph, replace = FALSE)
      cat(paste("      Nodos semilla principales (random):", paste(seed_nodes_to_test_as_primary, collapse=", "), "\n"))
    }
    
    if (length(seed_nodes_to_test_as_primary) == 0) {
      cat("      No se seleccionaron nodos semilla principales para este grafo. Saltando.\n")
      next
    }
    
    cl <- makeCluster(detectCores() - 1, type = "FORK") 
    registerDoParallel(cl)
    
    list_of_dfs_from_parallel_seeds <- foreach(
      current_primary_seed_id = seed_nodes_to_test_as_primary, 
      .combine = 'list', 
      .packages = c('igraph', 'dplyr'), 
      .export = c('sweep_homoph_parameter', 'get_complex_plot', 
                  'N_nodes_global', 'current_graph_obj_sim', 
                  'node_individual_thresholds_tau_frac_for_sim', # τ fraccional para get_complex_plot
                  'node_thresholds_count_for_plci_and_cluster', # τ conteo para clúster inicial
                  'node_mur_q_for_sim', 'node_degrees_for_sim',
                  'alpha_values_sim', 'homoph_values_sim'),
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
        all_innovation_iul_Gamma_values = alpha_values_sim, 
        all_social_distance_h_values = homoph_values_sim,   
        initial_infectors_vector_arg = initial_infectors_for_this_run # El clúster inicial
      )
      
      df_from_worker$graph_type <- current_graph_type_label
      df_from_worker$base_threshold_tau_frac <- current_threshold_base_tau_fractional 
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

cat("Todas las simulaciones completadas.\n")

# -----------------------------------------------------------------------------
# 4. Visualización Rápida de Resultados 
# -----------------------------------------------------------------------------
if (length(all_simulation_results_collection) > 0) {
  plots_per_tau_threshold <- list()
  
  # Definir nombres para la faceta si vas a usar múltiples tipos de grafos en un solo plot
  # graph_type_labels_for_plot <- c("Erdos-Renyi", "Scale-Free", "Small-World", "Small-World-SDA") # Ejemplo
  # names(graph_type_labels_for_plot) <- unique(sapply(all_simulation_results_collection, function(df) df$graph_type[1]))
  
  
  for (tau_label_for_plot_str in names(all_simulation_results_collection)) {
    data_for_this_tau_plot <- all_simulation_results_collection[[tau_label_for_plot_str]]
    
    # Asegurar que las columnas de resultados de get_complex_plot tengan los nombres correctos
    # (innovation_iul_Gamma, social_distance_h, seed, num_adopters)
    if (!is.null(data_for_this_tau_plot) && nrow(data_for_this_tau_plot) > 0 &&
        all(c("innovation_iul_Gamma", "num_adopters", "seed", "social_distance_h") %in% names(data_for_this_tau_plot))) {
      
      plot_title_text <- paste(current_graph_type_label, "- Umbral Base τ (frac):", tau_label_for_plot_str)
      
      data_for_this_tau_plot$social_distance_h_factor <- as.factor(sprintf("%.2f", data_for_this_tau_plot$social_distance_h))
      h_labels_for_legend <- sprintf("%.2f", unique(sort(data_for_this_tau_plot$social_distance_h)))
      
      fig_single_tau <- ggplot(data = data_for_this_tau_plot, 
                               aes(x = innovation_iul_Gamma, y = num_adopters / N_nodes_global, 
                                   color = social_distance_h_factor, group = social_distance_h_factor)) +
        # Quitar líneas individuales para claridad si hay muchas corridas
        # geom_line(aes(group = interaction(seed, graph_instance_idx)), linewidth = 0.3, alpha = 0.15) +
        stat_summary(fun = mean, geom = "line", linewidth = 1.0) +
        stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, aes(fill=social_distance_h_factor), show.legend = FALSE) + # Intervalo de confianza (opcional)
        scale_color_viridis_d(option = "plasma", name = "Distancia Social (h)", labels = h_labels_for_legend) +
        scale_fill_viridis_d(option = "plasma") + # Para el ribbon
        labs(x = expression(paste("Utilidad Intrínseca Innovación (", Gamma, ")")), 
             y = "Proporción Media de Adoptadores", 
             title = plot_title_text) +
        ylim(0,1) + 
        theme_minimal(base_size = 10) +
        theme(legend.position = "bottom", legend.direction = "horizontal", 
              plot.title = element_text(hjust = 0.5, size=11),
              strip.text = element_text(size=10)) # Para facet_wrap si se usa
      
      # Si quieres separar por graph_instance_idx (las 5 redes de la misma topología)
      # fig_single_tau <- fig_single_tau + facet_wrap(~graph_instance_idx, ncol=5)
      
      plots_per_tau_threshold[[plot_title_text]] <- fig_single_tau
      #print(fig_single_tau) 
      
    } else {
      cat(paste("No hay datos suficientes o faltan columnas para graficar el umbral τ:", tau_label_for_plot_str, "\n"))
      if(!is.null(data_for_this_tau_plot)) print(head(names(data_for_this_tau_plot))) # Ayuda a debuggear nombres de columnas
    }
  }
  
  if (length(plots_per_tau_threshold) > 1) {
    final_combined_plot_for_type <- wrap_plots(plots_per_tau_threshold, ncol = 1, guides = "collect") &
      theme(legend.position = "bottom")
    print(final_combined_plot_for_type)
    # ggsave(paste0("quick_plot_summary_", current_graph_type_label, ".pdf"), final_combined_plot_for_type, width=7, height=2.5*length(plots_per_tau_threshold), units="in")
  } else if (length(plots_per_tau_threshold) == 1) {
    # ggsave(paste0("quick_plot_summary_", current_graph_type_label, "_tau", names(plots_per_tau_threshold)[1] , ".pdf"), plots_per_tau_threshold[[1]], width=7, height=3, units="in")
  }
  
} else {
  cat("No se generaron resultados consolidados para visualizar.\n")
}
