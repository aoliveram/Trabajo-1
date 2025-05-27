# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script 2: simulation_functions.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Descripción: Contiene las funciones utilizadas por 'simulation_core.R'
#              para generar redes y ejecutar simulaciones de contagio.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(igraph)
library(doParallel) # para paralelización en sweep_homoph_parameter

# -----------------------------------------------------------------------------
# Generación de Redes
# -----------------------------------------------------------------------------

# Redes Small-World (Watts-Strogatz)
generate_small_world_networks <- function(N, num_edges, p_rewire) {
  small_world_graphs <- list()
  for (i in 1:5) { # Genera 5 instancias
    set.seed(i)
    k_neighbors <- round(2 * num_edges / N) 
    # Asegurar k_neighbors par y al menos 2
    if (k_neighbors < 2) k_neighbors <- 2
    if (k_neighbors %% 2 != 0) k_neighbors <- k_neighbors + 1
    
    graph <- sample_smallworld(dim = 1, size = N, nei = k_neighbors / 2, p = p_rewire)
    V(graph)$alpha <- runif(N, 0, 1) # Atributo 'alpha' para selectividad/homofilia
    small_world_graphs[[i]] <- graph
  }
  return(small_world_graphs)
}

# Redes Libres de Escala (Barabási-Albert)
generate_scale_free_networks <- function(N, m_param) {
  scale_free_graphs <- list()
  for (i in 1:5) {
    set.seed(i)
    graph <- barabasi.game(N, m = m_param, directed = FALSE)
    V(graph)$alpha <- runif(N, 0, 1)
    scale_free_graphs[[i]] <- graph
  }
  return(scale_free_graphs)
}

# Redes Erdös-Rényi G(N,M)
generate_erdos_renyi_networks <- function(N, num_edges_total) {
  erdos_renyi_graphs <- list()
  for (i in 1:5) {
    set.seed(i)
    graph <- erdos.renyi.game(N, num_edges_total, type = "gnm")
    V(graph)$alpha <- runif(N, 0, 1)
    erdos_renyi_graphs[[i]] <- graph
  }
  return(erdos_renyi_graphs)
}

# -----------------------------------------------------------------------------
# Funciones Auxiliares y de Simulación
# -----------------------------------------------------------------------------

State <- function(value = "State1") {
  if (!value %in% c("State1", "State2")) {
    stop("Estado no válido. Use 'State1' o 'State2'.")
  }
  return(value)
}

clustered_seeding <- function(current_seeds, graph_obj, num_total_seeds_needed) {
  # Asegurar que current_seeds sea numérico si viene de V(g)
  if (inherits(current_seeds, "igraph.vs")) {
    current_seeds <- as.numeric(current_seeds)
  }
  
  num_additional_seeds_needed <- num_total_seeds_needed - length(current_seeds)
  
  if (num_additional_seeds_needed > 0) {
    potential_new_seeds <- setdiff(as.numeric(V(graph_obj)), current_seeds)
    if (length(potential_new_seeds) < num_additional_seeds_needed) {
      # No hay suficientes nodos únicos para agregar, tomar todos los disponibles
      new_seeds_selected <- potential_new_seeds
      warning("No hay suficientes nodos únicos para agregar como semillas adicionales.")
    } else {
      new_seeds_selected <- sample(potential_new_seeds, num_additional_seeds_needed)
    }
    return(unique(c(current_seeds, new_seeds_selected)))
  } else {
    if (length(current_seeds) > num_total_seeds_needed && num_total_seeds_needed > 0) {
      
      # Para una lógica más general:
      # return(sample(current_seeds, num_total_seeds_needed))
      return(current_seeds) # Devolver las semillas actuales si ya son suficientes/más
    }
    return(unique(current_seeds))
  }
}
