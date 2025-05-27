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
                             node_individual_thresholds_tau, # τ_i ('thresholds_sim')
                             node_mur_q,                     # q_i (V(g)$alpha)
                             innovation_iul_Gamma,           # Γ ('target_alpha_1')
                             social_distance_h,              # h ('homoph_sensitivity')
                             initial_seed_nodes_vector) {    # Vector de nodos semilla iniciales
  
  # seed_node: el ID del nodo semilla principal (para tracking/output, no se usa directamente en la lógica de activación si ya está en initial_seed_nodes_vector)
  # N_nodes: número total de nodos
  # graph_obj: el objeto igraph
  # node_individual_thresholds_tau: vector de umbrales individuales τ_i para cada nodo
  # node_mur_q: vector de requerimientos de utilidad mínima q_i para cada nodo (V(g)$alpha)
  # innovation_iul_Gamma: escalar, utilidad intrínseca de la innovación Γ
  # social_distance_h: escalar, umbral de cercanía social h
  # initial_seed_nodes_vector: vector con los índices de los nodos que inician activados
  
  node_states_activated <- rep("State1", N_nodes)
  if (length(initial_seed_nodes_vector) > 0) {
    node_states_activated[initial_seed_nodes_vector] <- "State2"
  } else {
    # Si no hay semillas iniciales, no puede haber propagación.
    # Devolver 0 adoptadores o manejar como error.
    # Por ahora, asumimos que siempre habrá al menos una semilla.
    # Si initial_seed_nodes_vector es vacío, la función devolverá num_final_adopters = 0
  }
  
  adj_matrix <- as_adjacency_matrix(graph_obj, sparse = FALSE) 
  
  spread_continues <- TRUE
  simulation_step <- 0 # Empezar en 0, primer chequeo es paso 1
  
  # Bucle principal de simulación (mientras haya cambios)
  while (spread_continues) {
    simulation_step <- simulation_step + 1
    newly_infected_this_step_idx <- c()
    
    # Nodos actualmente activos
    current_activated_nodes_idx <- which(node_states_activated == "State2")
    if (length(current_activated_nodes_idx) == 0 && simulation_step > 1) { 
      # Si no hay nadie activo (y no es el inicio), no hay más propagación
      spread_continues <- FALSE
      next
    }
    
    # Iterar sobre todos los nodos INACTIVOS
    inactive_nodes_idx <- which(node_states_activated == "State1")
    
    if (length(inactive_nodes_idx) == 0) { # Todos están activos
      spread_continues <- FALSE
      next
    }
    
    for (i in inactive_nodes_idx) {
      # El nodo 'i' está inactivo.
      # ¿Se entera de la innovación? (i.e., ¿tiene al menos un vecino activo?)
      neighbors_of_i_idx <- which(adj_matrix[i, ] > 0)
      active_neighbors_of_i_idx <- intersect(neighbors_of_i_idx, current_activated_nodes_idx)
      
      if (length(active_neighbors_of_i_idx) > 0) { # El nodo 'i' se entera
        
        # Aplicar Regla 3: Elección Racional
        if (node_mur_q[i] <= innovation_iul_Gamma) {
          newly_infected_this_step_idx <- c(newly_infected_this_step_idx, i)
          # Pasa a 'continue' para no evaluar la regla 4 para este nodo en este paso
          next 
        }
        
        # Aplicar Regla 4: Influencia Social (si q_i > Γ)
        
        # Calcular Exposición Homofílica Ẽ_i (vecinos ACTIVOS y SOCIALMENTE CERCANOS)
        exposure_Ei_tilde <- 0
        
        if (length(active_neighbors_of_i_idx) > 0) {
          for (j in active_neighbors_of_i_idx) {
            # Calcular distancia social d_ij (basada en la diferencia de q_i, que es node_mur_q)
            social_dist_ij <- abs(node_mur_q[i] - node_mur_q[j])
            
            if (social_dist_ij <= social_distance_h) {
              # Ẽ_i es la *cuenta* de vecinos activos y socialmente cercanos.
              exposure_Ei_tilde <- exposure_Ei_tilde + 1
            }
          }
        }
        
        # Comparar con el umbral individual τ_i
        if (exposure_Ei_tilde >= node_individual_thresholds_tau[i]) {
          newly_infected_this_step_idx <- c(newly_infected_this_step_idx, i)
        }
      } # Fin if (el nodo 'i' se entera)
    } # Fin for (iterar sobre nodos inactivos)
    
    # Actualizar estados y verificar si hubo cambios
    if (length(newly_infected_this_step_idx) > 0) {
      node_states_activated[unique(newly_infected_this_step_idx)] <- "State2"
      # spread_continues sigue TRUE
    } else {
      spread_continues <- FALSE # No hubo nuevos infectados en este paso
    }
    
    if (all(node_states_activated == "State2")) { # Todos infectados
      spread_continues <- FALSE
    }
    
    if (simulation_step > (N_nodes + 5) ) { # Salvaguarda
      #warning(paste("Simulación para semilla principal", seed_node, "detenida por exceso de pasos."))
      spread_continues <- FALSE
    }
  } # Fin while (spread_continues)
  
  num_final_adopters <- sum(node_states_activated == "State2")
  
  # Devolver los resultados para esta corrida de simulación
  # 'seed' aquí se refiere al parámetro 'seed_node' que identifica la corrida, no necesariamente el único nodo semilla.
  return(list(
    alpha_1 = innovation_iul_Gamma, 
    homophily = social_distance_h, 
    seed = seed_node, # El ID de la semilla principal de esta simulación particular
    num_adopters = num_final_adopters, 
    num_steps = simulation_step
  ))
}

