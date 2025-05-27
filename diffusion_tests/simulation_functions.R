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

# En simulation_functions.R

get_complex_plot <- function(seed_node, N_nodes, graph_obj, 
                             node_individual_thresholds_tau, # τ_i 
                             node_mur_q,                     # q_i (V(g)$alpha)
                             innovation_iul_Gamma,           # Γ 
                             social_distance_h,              # h 
                             initial_seed_nodes_vector) {    # Vector de nodos semilla iniciales
  
  node_states_activated <- rep("State1", N_nodes)
  if (length(initial_seed_nodes_vector) > 0) {
    node_states_activated[initial_seed_nodes_vector] <- "State2"
  }
  
  adj_matrix <- as_adjacency_matrix(graph_obj, sparse = FALSE) 
  
  spread_continues <- TRUE
  simulation_step <- 0
  
  while (spread_continues) {
    simulation_step <- simulation_step + 1
    newly_infected_this_step_idx <- c()
    current_activated_nodes_idx <- which(node_states_activated == "State2")
    
    if (length(current_activated_nodes_idx) == 0 && simulation_step > 1) {
      spread_continues <- FALSE
      next
    }
    
    inactive_nodes_idx <- which(node_states_activated == "State1")
    if (length(inactive_nodes_idx) == 0) {
      spread_continues <- FALSE
      next
    }
    
    for (i in inactive_nodes_idx) {
      neighbors_of_i_idx <- which(adj_matrix[i, ] > 0)
      active_neighbors_of_i_idx <- intersect(neighbors_of_i_idx, current_activated_nodes_idx)
      
      if (length(active_neighbors_of_i_idx) > 0) { # El nodo 'i' se entera
        if (node_mur_q[i] <= innovation_iul_Gamma) { # Regla 3: Elección Racional
          newly_infected_this_step_idx <- c(newly_infected_this_step_idx, i)
          next 
        }
        
        # Regla 4: Influencia Social (si q_i > Γ)
        # Calcular Exposición Homofílica Ẽ_i según Ecuación (4)
        
        # Numerador: Suma de (X_ij * a_j) donde a_j implica activo Y d_ij <= h
        numerator_Ei_tilde <- 0
        # Denominador: Suma de X_ij (es decir, el grado del nodo i)
        denominator_Ei_tilde <- length(neighbors_of_i_idx) # Grado del nodo i
        
        if (denominator_Ei_tilde > 0 && length(active_neighbors_of_i_idx) > 0) {
          for (j in active_neighbors_of_i_idx) { # Iterar solo sobre vecinos ACTIVOS
            social_dist_ij <- abs(node_mur_q[i] - node_mur_q[j]) # d_ij
            if (social_dist_ij <= social_distance_h) {
              # El vecino 'j' es activo Y socialmente cercano.
              # X_ij es 1 (porque es vecino), a_j (tilde) es 1 (porque es activo y cercano)
              numerator_Ei_tilde <- numerator_Ei_tilde + 1 
            }
          }
          
          exposure_Ei_tilde_value <- numerator_Ei_tilde / denominator_Ei_tilde
          
          # Comparar con el umbral individual τ_i
          # La tesis dice τ_i ≤ Ẽ_i. En tu código, los umbrales son generalmente enteros (conteo).
          # Si Ẽ_i es una proporción (0 a 1) y τ_i es un conteo, la comparación directa no es correcta.
          # Si τ_i también es una proporción (0 a 1), entonces sí.
          # Tu código original y `threshold_values_list_sim` sugieren que τ_i son fracciones del grado,
          # o valores pequeños como 0.10, 0.15. Si son fracciones (0 a 1), la comparación es válida.
          if (exposure_Ei_tilde_value >= node_individual_thresholds_tau[i]) {
            newly_infected_this_step_idx <- c(newly_infected_this_step_idx, i)
          }
        } else {
          # Si no tiene vecinos (denominator_Ei_tilde == 0) o no tiene vecinos activos,
          # exposure_Ei_tilde_value será 0 o NaN, y no se activará por esta regla.
        }
      } 
    } 
    
    if (length(newly_infected_this_step_idx) > 0) {
      node_states_activated[unique(newly_infected_this_step_idx)] <- "State2"
    } else {
      spread_continues <- FALSE 
    }
    
    if (all(node_states_activated == "State2")) {
      spread_continues <- FALSE
    }
    if (simulation_step > (N_nodes + 5) ) {
      spread_continues <- FALSE
    }
  } 
  
  num_final_adopters <- sum(node_states_activated == "State2")
  return(list(
    innovation_iul_Gamma = innovation_iul_Gamma, # Cambiado de alpha_1
    social_distance_h = social_distance_h,    # Cambiado de homophily
    seed = seed_node, 
    num_adopters = num_final_adopters, 
    num_steps = simulation_step
  ))
}
