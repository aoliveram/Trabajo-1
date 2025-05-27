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
          
          # Comparar con el umbral individual τ_i (fraccionario)
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

sweep_homoph_parameter <- function(primary_seed_id_arg, 
                                   N_nodes_arg, 
                                   graph_obj_arg, 
                                   node_individual_thresholds_tau_arg, # Vector τ_i
                                   node_mur_q_arg,                     # Vector q_i
                                   all_innovation_iul_Gamma_values,  # Vector de Γ a probar
                                   all_social_distance_h_values,     # Vector de h a probar
                                   initial_infectors_vector_arg) {   # Semillas para ESTA simulación
  
  results_for_this_primary_seed <- list()
  
  for (current_Gamma in all_innovation_iul_Gamma_values) {
    for (current_h in all_social_distance_h_values) {
      
      sim_output <- get_complex_plot(
        seed_node = primary_seed_id_arg, # Para identificar la corrida
        N_nodes = N_nodes_arg,
        graph_obj = graph_obj_arg,
        node_individual_thresholds_tau = node_individual_thresholds_tau_arg,
        node_mur_q = node_mur_q_arg,
        innovation_iul_Gamma = current_Gamma,
        social_distance_h = current_h,
        initial_seed_nodes_vector = initial_infectors_vector_arg
      )
      results_for_this_primary_seed[[length(results_for_this_primary_seed) + 1]] <- as.data.frame(sim_output)
    }
  }
  
  df_for_this_primary_seed <- bind_rows(results_for_this_primary_seed)
  return(df_for_this_primary_seed)
}

# Funciones de Centralidad (Mantenidas con cambio de nombre)
get_simple_centralities <- function(graph_obj) {
  centrality_df <- data.frame(
    seed = as.numeric(V(graph_obj)), 
    degree = as.numeric(degree(graph_obj)), 
    betweenness = as.numeric(betweenness(graph_obj)), 
    eigen = as.numeric(eigen_centrality(graph_obj)$vector)
  ) #centrality_df$percolation <- NA
  return(centrality_df)
}

get_PLci <- function(seed_node, N_nodes, graph_obj, adj_matrix, node_thresholds, num_initial_cluster_seeds_for_node, model_output_list_placeholder) {
  
  adj_matrix_simulation_run_PLci <- matrix(0, nrow = N_nodes, ncol = N_nodes)
  
  # Selección de semillas para esta simulación de PLci (seed + sus num_initial_cluster_seeds_for_node-1 vecinos/random)
  # num_initial_cluster_seeds_for_node es el umbral para ESE nodo.
  initial_plci_seeds <- c(seed_node) # Empezar con la semilla
  if(num_initial_cluster_seeds_for_node > 1) { # Si se necesitan más semillas en el cluster
    neighbor_indices <- as.numeric(neighbors(graph_obj, seed_node, mode="total"))
    if (length(neighbor_indices) >= (num_initial_cluster_seeds_for_node -1) ) {
      initial_plci_seeds <- c(initial_plci_seeds, sample(neighbor_indices, num_initial_cluster_seeds_for_node -1))
    } else {
      initial_plci_seeds <- c(initial_plci_seeds, neighbor_indices)
      # Si aún faltan, agregar aleatorios (excluyendo los ya seleccionados)
      needed_more <- num_initial_cluster_seeds_for_node - length(initial_plci_seeds)
      if (needed_more > 0) {
        potential_others <- setdiff(1:N_nodes, initial_plci_seeds)
        if (length(potential_others) >= needed_more) {
          initial_plci_seeds <- c(initial_plci_seeds, sample(potential_others, needed_more))
        } else {
          initial_plci_seeds <- c(initial_plci_seeds, potential_others)
        }
      }
    }
  }
  initial_plci_seeds <- unique(initial_plci_seeds)
  
  activated_plci <- rep(FALSE, N_nodes)
  if(length(initial_plci_seeds) > 0) {
    activated_plci[initial_plci_seeds] <- TRUE
    adj_matrix_simulation_run_PLci[initial_plci_seeds, ] <- 1
    adj_matrix_simulation_run_PLci[, initial_plci_seeds] <- 1
  }
  current_adj_matrix_active_links_PLci <- adj_matrix * adj_matrix_simulation_run_PLci # Aristas activas iniciales
  
  spread_continues_plci <- TRUE
  num_steps_plci <- 0
  while (spread_continues_plci) {
    num_steps_plci <- num_steps_plci + 1
    influence_plci <- colSums(current_adj_matrix_active_links_PLci)
    nodes_passing_threshold_plci <- influence_plci >= node_thresholds
    
    newly_activated_this_step_plci <- which(nodes_passing_threshold_plci & !activated_plci)
    
    if (length(newly_activated_this_step_plci) == 0) {
      spread_continues_plci <- FALSE
    } else {
      activated_plci[newly_activated_this_step_plci] <- TRUE
      # Actualizar adj_matrix_simulation_run_PLci para las nuevas aristas activas
      # Solo las aristas de los *nuevos* adoptadores a todos los demás, y viceversa
      # O más simple: reconstruir desde todos los activos
      all_activated_indices_plci <- which(activated_plci)
      adj_matrix_simulation_run_PLci[,] <- 0 # Reset
      if(length(all_activated_indices_plci)>0){
        adj_matrix_simulation_run_PLci[all_activated_indices_plci, ] <- 1
        adj_matrix_simulation_run_PLci[, all_activated_indices_plci] <- 1
      }
      current_adj_matrix_active_links_PLci <- adj_matrix * adj_matrix_simulation_run_PLci
    }
    if(all(activated_plci) || num_steps_plci > N_nodes + 5) spread_continues_plci <- FALSE
  }
  
  num_final_adopters_plci <- sum(activated_plci)
  
  complex_spread_graph <- graph_from_adjacency_matrix(current_adj_matrix_active_links_PLci, mode = "undirected")
  
  reachable_nodes <- which(activated_plci)
  if (seed_node %in% reachable_nodes && length(reachable_nodes) > 0) {
    distances_from_seed <- distances(complex_spread_graph, v = seed_node, to = reachable_nodes, mode = "out")
    # Reemplazar Inf (no alcanzable dentro del cluster de activados, no debería pasar si seed está) por 0 o N
    # distances_from_seed[is.infinite(distances_from_seed)] <- N_nodes 
    # Si solo queremos la media de los alcanzados:
    distances_from_seed[is.infinite(distances_from_seed)] <- NA
    
    
    PLci_value <- mean(distances_from_seed, na.rm = TRUE)
    if(is.nan(PLci_value)) PLci_value <- N_nodes # Si no se alcanzó a nadie más que la semilla o solo la semilla y no hay distancias
  } else {
    PLci_value <- N_nodes # Penalización si la semilla no alcanza a nadie o no se activa
  }
  
  # Devolver lista estructurada como en tu código original de get_complex
  return(c(seed_node, N_nodes, num_initial_cluster_seeds_for_node, num_final_adopters_plci, PLci_value))
}
