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
