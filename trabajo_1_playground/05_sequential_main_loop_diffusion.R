# -----------------------------------------------------------------------------
# 3. Bucle Principal de Simulación (VERSIÓN SECUENCIAL)
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
    
    node_mur_q_for_sim <- V(current_graph_obj_sim)$q_i 
    node_degrees_for_sim <- igraph::degree(current_graph_obj_sim)
    
    # Calcular la matriz de disimilitud de Gower
    if (ATP_NETWORK || all(c("age", "educ_num", "race", "relig", "sex") %in% igraph::vertex_attr_names(current_graph_obj_sim))) {
      attributes_for_distance <- data.frame(
        age = V(current_graph_obj_sim)$age,
        educ_num = V(current_graph_obj_sim)$educ_num,
        race = as.factor(V(current_graph_obj_sim)$race),
        relig = as.factor(V(current_graph_obj_sim)$relig),
        sex = as.factor(V(current_graph_obj_sim)$sex)
      )
      social_distance_matrix_d_ij <- as.matrix(daisy(attributes_for_distance, metric = "gower"))
    } else {
      cat("      ADVERTENCIA: Atributos para Gower no encontrados.\n")
      break
    }
    
    # Calcular umbrales individuales τ_i
    node_individual_thresholds_tau_frac_for_sim <- rep(current_threshold_base_tau_fractional, N_nodes_global)
    node_individual_thresholds_tau_frac_for_sim[node_individual_thresholds_tau_frac_for_sim <= 0 & current_threshold_base_tau_fractional > 0] <- 1e-6
    node_individual_thresholds_tau_frac_for_sim[node_individual_thresholds_tau_frac_for_sim > 1] <- 1.0
    
    node_thresholds_count_for_plci_and_cluster <- round(node_individual_thresholds_tau_frac_for_sim * node_degrees_for_sim)
    node_thresholds_count_for_plci_and_cluster[node_thresholds_count_for_plci_and_cluster <= 0 & node_degrees_for_sim > 0] <- 1 # Mínimo 1 si el grado es > 0
    node_thresholds_count_for_plci_and_cluster[node_degrees_for_sim == 0] <- 0 # Para nodos aislados, el umbral de conteo debe ser 0
    
    seed_nodes_to_test_as_primary <- c()
    if (SEEDING_STRATEGY == "PLci_top") {
      adj_mat_for_plci_calc <- as.matrix(as_adjacency_matrix(current_graph_obj_sim))
      # Asegurarse de que N_nodes_global coincida con el grafo actual si los grafos pueden variar en tamaño
      current_N_nodes <- vcount(current_graph_obj_sim)
      
      plci_values <- sapply(1:current_N_nodes, function(node_idx_for_plci) {
        # Verificar que el tamaño del clúster inicial para PLCI no sea problemático
        num_initial_seeds_for_node_plci <- node_thresholds_count_for_plci_and_cluster[node_idx_for_plci]
        if (num_initial_seeds_for_node_plci > current_N_nodes) num_initial_seeds_for_node_plci <- current_N_nodes
        if (num_initial_seeds_for_node_plci < 0 ) num_initial_seeds_for_node_plci <- 0 # O 1 si se espera al menos 1

        get_PLci(
          seed_node = node_idx_for_plci, 
          N_nodes = current_N_nodes, # Usar N_nodes del grafo actual
          graph_obj = current_graph_obj_sim, 
          adj_matrix = adj_mat_for_plci_calc, 
          node_thresholds = node_thresholds_count_for_plci_and_cluster, 
          num_initial_cluster_seeds_for_node = node_thresholds_count_for_plci_and_cluster[node_idx_for_plci],
          model_output_list_placeholder = NULL 
        )[5] 
      })
      seed_nodes_to_test_as_primary <- order(as.numeric(plci_values), decreasing = TRUE)[1:min(current_N_nodes, NUM_INITIAL_SEEDS)]
      cat(paste("      Nodos semilla principales (PLci_top):", paste(seed_nodes_to_test_as_primary, collapse=", "), "\n"))
    } else if (SEEDING_STRATEGY == "random") {
      set.seed(graph_idx * sum(as.numeric(charToRaw(T_str_sim))) * 7) 
      current_N_nodes <- vcount(current_graph_obj_sim)
      if (current_N_nodes > 0) {
        seed_nodes_to_test_as_primary <- as.numeric(sample(V(current_graph_obj_sim), min(NUM_INITIAL_SEEDS, current_N_nodes), replace = FALSE))
      } else {
        seed_nodes_to_test_as_primary <- c()
      }
      cat(paste("      Nodos semilla principales (random):", paste(seed_nodes_to_test_as_primary, collapse=", "), "\n"))
    }
    
    if (length(seed_nodes_to_test_as_primary) == 0) {
      cat("      No se seleccionaron nodos semilla principales para este grafo. Saltando.\n")
      next # Saltar al siguiente grafo si no hay semillas
    }
    
    # Lista para almacenar los data.frames de resultados para esta instancia de grafo
    list_of_dfs_from_sequential_seeds <- list()
    
    # Bucle secuencial sobre los nodos semilla a probar
    for (current_primary_seed_id in seed_nodes_to_test_as_primary) {
      
      cat(paste("        Procesando semilla principal:", current_primary_seed_id, "\n")) # Log
      
      current_N_nodes <- vcount(current_graph_obj_sim) # Reafirmar N_nodes para el grafo actual
      
      num_seeds_for_initial_cluster = node_thresholds_count_for_plci_and_cluster[current_primary_seed_id]
      num_seeds_for_initial_cluster = min(num_seeds_for_initial_cluster, current_N_nodes, node_degrees_for_sim[current_primary_seed_id] + 1)
      if(num_seeds_for_initial_cluster < 1 && current_N_nodes > 0) num_seeds_for_initial_cluster <- 1
      if(current_N_nodes == 0) num_seeds_for_initial_cluster <- 0
      
      
      initial_infectors_for_this_run <- c()
      if (current_N_nodes > 0) { # Solo intentar seleccionar semillas si hay nodos
        initial_infectors_for_this_run <- c(current_primary_seed_id)
        
        if (num_seeds_for_initial_cluster > 1) {
          neighbors_of_primary <- as.numeric(neighbors(current_graph_obj_sim, current_primary_seed_id, mode="total"))
          num_additional_needed = num_seeds_for_initial_cluster - 1
          
          if (length(neighbors_of_primary) >= num_additional_needed) {
            initial_infectors_for_this_run <- c(initial_infectors_for_this_run, sample(neighbors_of_primary, num_additional_needed))
          } else { 
            initial_infectors_for_this_run <- c(initial_infectors_for_this_run, neighbors_of_primary)
            still_needed_more <- num_seeds_for_initial_cluster - length(initial_infectors_for_this_run)
            if (still_needed_more > 0) {
              potential_others <- setdiff(1:current_N_nodes, initial_infectors_for_this_run)
              if (length(potential_others) > 0) {
                initial_infectors_for_this_run <- c(initial_infectors_for_this_run, sample(potential_others, min(length(potential_others), still_needed_more)))
              }
            }
          }
        }
        initial_infectors_for_this_run <- unique(initial_infectors_for_this_run)
      }
      
      # Capturar errores potenciales de sweep_homoph_parameter
      df_from_worker <- tryCatch({
        temp_df <- sweep_homoph_parameter(
          primary_seed_id_arg = current_primary_seed_id,
          N_nodes_arg = current_N_nodes, # Usar N_nodes del grafo actual
          graph_obj_arg = current_graph_obj_sim,
          node_individual_thresholds_tau_arg = node_individual_thresholds_tau_frac_for_sim,
          node_mur_q_arg = node_mur_q_for_sim,
          all_innovation_iul_Gamma_values = mur_values_sim, 
          all_social_distance_h_values = homoph_values_sim,   
          initial_infectors_vector_arg = initial_infectors_for_this_run,
          d_ij_matrix = social_distance_matrix_d_ij
        )
        
        temp_df$graph_type <- current_graph_type_label
        temp_df$base_threshold_tau_frac <- current_threshold_base_tau_fractional 
        temp_df$graph_instance_idx <- graph_idx
        temp_df # Devolver el data.frame
      }, error = function(e) {
        cat(paste("        ERROR al procesar semilla", current_primary_seed_id, ":", e$message, "\n"))
        return(NULL) # Devolver NULL o un objeto de error si se prefiere
      })
      
      # Añadir el resultado (data.frame o NULL) a la lista
      if (!is.null(df_from_worker)) {
        list_of_dfs_from_sequential_seeds[[length(list_of_dfs_from_sequential_seeds) + 1]] <- df_from_worker
      }
    } # Fin del bucle secuencial sobre seed_nodes_to_test_as_primary
    
    
    # Filtrar y combinar resultados (similar a la lógica post-paralelización)
    final_dfs_to_combine_for_this_graph_instance <- list()
    successful_runs_count_for_filter = 0
    
    for(df_from_a_seed_run in list_of_dfs_from_sequential_seeds){
      # Ya hemos manejado errores con tryCatch, así que df_from_a_seed_run será un data.frame o NULL
      if (is.null(df_from_a_seed_run) || !is.data.frame(df_from_a_seed_run) || nrow(df_from_a_seed_run) == 0) {
        next
      }
      
      process_this_df = FALSE
      # N_nodes_global aquí se refiere al tamaño original esperado,
      # pero num_adopters viene de la simulación en current_graph_obj_sim.
      # Para la proporción, usar el tamaño del grafo actual (current_N_nodes).
      current_N_nodes_for_proportion <- vcount(current_graph_obj_sim) 
      
      if (SEEDING_STRATEGY == "PLci_top") {
        if (successful_runs_count_for_filter < NUM_SUCCESSFUL_SIMS_PER_GRAPH) {
          if (any(df_from_a_seed_run$num_adopters / current_N_nodes_for_proportion > num_adopters_min_sim, na.rm=TRUE)) {
            process_this_df = TRUE
          }
        }
      } else { # Para "random"
        if (successful_runs_count_for_filter < NUM_SUCCESSFUL_SIMS_PER_GRAPH) {
          process_this_df = TRUE
        }
      }
      
      if(process_this_df){
        final_dfs_to_combine_for_this_graph_instance[[length(final_dfs_to_combine_for_this_graph_instance) + 1]] <- df_from_a_seed_run
        successful_runs_count_for_filter <- successful_runs_count_for_filter + 1
      }
    }
    
    if (length(final_dfs_to_combine_for_this_graph_instance) > 0) {
      results_for_this_base_tau_and_all_graphs[[paste0("graph_", graph_idx)]] <- bind_rows(final_dfs_to_combine_for_this_graph_instance)
    } else {
      cat(paste("    No se obtuvieron resultados válidos para el grafo #", graph_idx, "y umbral τ", T_str_sim, "\n"))
    }
  } # Fin del bucle sobre graph_idx
  
  if (length(results_for_this_base_tau_and_all_graphs) > 0) {
    all_simulation_results_collection[[T_str_sim]] <- bind_rows(results_for_this_base_tau_and_all_graphs)
  } else {
    cat(paste("  No se obtuvieron resultados para el umbral τ", T_str_sim, "\n"))
  }
} # Fin del bucle sobre current_threshold_base_tau_fractional 


# Guardamos
saveRDS(all_simulation_results_collection, "trabajo_1_files/diffusion_ATP_nets_2.rds")
