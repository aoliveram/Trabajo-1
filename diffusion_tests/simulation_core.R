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
  
  # graphs_list_to_simulate <- generate_scale_free_networks(N_nodes_sim, m_sim)
  # current_graph_type_label <- "Scale-Free"
  
  graphs_list_to_simulate <- generate_erdos_renyi_networks(N_nodes_sim, num_edges_sim)
  current_graph_type_label <- "Erdos-Renyi"
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
all_simulation_results <- list() # almacena los resultados

cat(paste("Iniciando simulaciones para tipo de red:", current_graph_type_label, "\n"))

for (current_threshold_base in threshold_values_list_sim) {
  T_str_sim <- as.character(current_threshold_base)
  cat(paste("  Corriendo para Umbral Base:", T_str_sim, "\n"))
  
  simulation_results_for_threshold <- list()
  
  for (graph_idx in seq_along(graphs_list_to_simulate)) {
    current_graph <- graphs_list_to_simulate[[graph_idx]]
    cat(paste("    Procesando grafo #", graph_idx, "de 5...\n"))
    
    # Calcular umbrales para los nodos
    thresholds_sim <- rep(current_threshold_base, N_nodes_global)
    if (T_type_sim == "frac") {
      degrees <- degree(current_graph, mode = "total")
      thresholds_sim <- round(thresholds_sim * degrees)
    }
    # Asegurar que los umbrales no sean cero o negativos si son absolutos
    thresholds_sim[thresholds_sim <= 0] <- 1 
    
    # num_seeds_to_add se usa en get_PLci para determinar el "tamaño" inicial del cluster de siembra
    num_seeds_to_add_sim <- thresholds_sim - 1
    num_seeds_to_add_sim[num_seeds_to_add_sim < 0] <- 0
    
    # Selección de nodos semilla
    seeds_to_run_simulation <- c()
    if (SEEDING_STRATEGY == "PLci_top") {
      # Calcular PLci para todos los nodos
      # Nota: get_PLci espera una lista 'model_output_list' que no se usa realmente.
      # Se puede simplificar en la función o pasar NULL si no causa error.
      temp_adj_matrix <- as.matrix(as_adjacency_matrix(current_graph))
      complex_centrality_values <- sapply(1:N_nodes_global, function(node_idx) {
        # Asegúrate de que 'model_output_list' no sea un problema aquí o modifica get_PLci
        get_PLci(node_idx, N_nodes_global, current_graph, temp_adj_matrix, thresholds_sim, num_seeds_to_add_sim, model_output_list = NULL)[5] # El 5to elemento es PLci
      })
      top_PLci_nodes_indices <- order(as.numeric(complex_centrality_values), decreasing = TRUE)[1:NUM_TOP_PLCI_NODES_TO_TEST]
      seeds_to_run_simulation <- top_PLci_nodes_indices
      cat(paste("      Nodos semilla (PLci_top):", paste(seeds_to_run_simulation, collapse=", "), "\n"))
      
    } else if (SEEDING_STRATEGY == "random") {
      # Para robustez con semillas aleatorias, idealmente se harían MÚLTIPLES corridas aquí
      # y se promediarían, o se tomaría una muestra representativa.
      # Por simplicidad para testeo rápido, solo tomamos un conjunto.
      set.seed(graph_idx * current_threshold_base * 100) # Semilla variable
      seeds_to_run_simulation <- sample(V(current_graph), NUM_SUCCESSFUL_SIMS_PER_GRAPH) # O un número fijo de nodos aleatorios a probar
      cat(paste("      Nodos semilla (random):", paste(seeds_to_run_simulation, collapse=", "), "\n"))
    }
    
    successful_sim_count_for_graph <- 0
    
    for (seed_node_sim in seeds_to_run_simulation) {
      if (successful_sim_count_for_graph >= NUM_SUCCESSFUL_SIMS_PER_GRAPH && SEEDING_STRATEGY == "PLci_top") {
        break # Ya tenemos suficientes simulaciones exitosas para esta estrategia
      }
      
      cat(paste("        Simulando con semilla:", seed_node_sim, "...\n"))
      
      # Ejecutar el barrido de parámetros para esta semilla
      # El parámetro 'plot' se quita de la llamada y de la función sweep_homoph_parameter
      # 'file_dir' y 'centrality' tampoco son necesarios para el testeo rápido sin guardado.
      # 'i' en la llamada original a sweep_homoph_parameter era el índice del grafo, ahora es graph_idx
      temp_results_df <- sweep_homoph_parameter(
        seed_node = seed_node_sim, 
        homoph_values = homoph_values_sim, 
        N = N_nodes_global, 
        graph = current_graph, 
        thresholds = thresholds_sim, 
        num_seeds_to_add = num_seeds_to_add_sim, 
        alpha_values = alpha_values_sim,
        graph_idx_for_saving = graph_idx # 'i' en tu código original
      )
      
      # Añadir info de la corrida actual
      temp_results_df$graph_type <- current_graph_type_label
      temp_results_df$base_threshold <- current_threshold_base
      temp_results_df$graph_instance_idx <- graph_idx
      
      # Verificar si la simulación fue "exitosa" (para estrategia PLci_top)
      # o simplemente guardar los resultados (para estrategia random)
      if (SEEDING_STRATEGY == "PLci_top") {
        if (any(temp_results_df$num_adopters / N_nodes_global > num_adopters_min_sim)) {
          simulation_results_for_threshold[[length(simulation_results_for_threshold) + 1]] <- temp_results_df
          successful_sim_count_for_graph <- successful_sim_count_for_graph + 1
          cat("          Simulación exitosa guardada.\n")
        } else {
          cat("          Simulación no exitosa (no superó el mínimo de adoptadores).\n")
        }
      } else { # Para random, guardamos todas las corridas
        simulation_results_for_threshold[[length(simulation_results_for_threshold) + 1]] <- temp_results_df
        successful_sim_count_for_graph <- successful_sim_count_for_graph + 1 # Contamos todas
      }
    } # Fin bucle seed_node_sim
  } # Fin bucle graph_idx
  
  # Combinar resultados para este umbral y añadirlos a la lista general
  if (length(simulation_results_for_threshold) > 0) {
    combined_df_for_threshold <- bind_rows(simulation_results_for_threshold)
    all_simulation_results[[T_str_sim]] <- combined_df_for_threshold
  } else {
    cat(paste("    No se obtuvieron resultados para el umbral:", T_str_sim, "\n"))
  }
  
} # Fin bucle current_threshold_base

cat("Todas las simulaciones completadas.\n")