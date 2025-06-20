# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Script to simulate complex contagion, identify phase transitions,
#              and visualize them as heatmaps of IUL vs. h.
#              Sweeps over different MEANS and STANDARD DEVIATIONS for the
#              normal threshold distribution.
#              Each simulation run uses a new network instance and a new random seed.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(igraph)
library(doParallel)
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(intergraph)
library(cluster)

# Load the script with simulation functions
source("diffusion_tests/simulation_functions.R") # Ensure this path is correct

# -----------------------------------------------------------------------------
# 1. Core Parameters & Setup
# -----------------------------------------------------------------------------
cat("Setting up core parameters...\n")

CURRENT_GRAPH_TYPE_LABEL <- "ATP-net"
NETWORKS_DIR <- "trabajo_1_files/ATP_network_ergm/"
# Ensure you have 100 network files if NUM_SEED_RUNS_TOTAL is 100
# Or adjust NUM_NETWORK_INSTANCES_AVAILABLE if you have fewer
NUM_NETWORK_INSTANCES_AVAILABLE <- 96 # Max number of unique network files
NUM_SEED_RUNS_TOTAL <- 96 # Total number of simulation runs (network+seed pairs) per (Mean, SD) combo
N_NODES_GLOBAL <- 1000

# Means for Normal Threshold Distribution (τ_i) - INNER SWEEP
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6)

# Standard Deviations for Normal Threshold Distribution (τ_i) - OUTER SWEEP
TAU_NORMAL_SD_SWEEP_LIST <- c(0.08, 0.12, 0.16, 0.20) # Your point 2

SEEDING_STRATEGY_FIXED <- "eigen" #"central" #"random"
PHASE_TRANSITION_THRESHOLD_JUMP <- 0.50
NUM_CORES_TO_USE <- 8

# Output directories
RESULTS_DIR <- "trabajo_1_files/diffusion_simulation_files/"
PLOTS_DIR <- "trabajo_1_plots/diffusion_simulation_plots/"
# dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
# dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. Parameters to Sweep (IUL and h - within each simulation run)
# -----------------------------------------------------------------------------
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP <- seq(0/12, 12/12, by = 1/12)

# -----------------------------------------------------------------------------
# 3. Main Simulation Loop (Outer loop for SD, Inner loop for Mean)
# -----------------------------------------------------------------------------
cat("Starting grand simulation sweep...\n")

# This will store all heatmap dataframes, one for each SD
all_sds_heatmap_data_list <- list()
# This will store all raw results, one list element per SD, 
# where each element is itself a list of results per Mean
all_sds_raw_results_list <- list()

total_sd_iterations <- length(TAU_NORMAL_SD_SWEEP_LIST)

time_init <- Sys.time()
for (sd_idx in 1:total_sd_iterations) {
  current_tau_sd <- TAU_NORMAL_SD_SWEEP_LIST[sd_idx]
  
  cat(paste0("\n####################################################################\n"))
  cat(paste0("Processing for Tau Distribution Standard Deviation (σ): ", current_tau_sd, "\n"))
  cat(paste0(" (SD Iteration ", sd_idx, " of ", total_sd_iterations, ")\n"))
  cat(paste0("####################################################################\n"))
  
  # This will store raw results for the current SD, one list element per Mean
  current_sd_all_means_raw_results <- list()
  # This will store heatmap data for the current SD, one list element per Mean (for faceting)
  current_sd_all_means_heatmap_data <- list()
  
  total_mean_iterations <- length(THRESHOLD_MEAN_SWEEP_LIST)
  for (mean_idx in 1:total_mean_iterations) {
    current_threshold_mean <- THRESHOLD_MEAN_SWEEP_LIST[mean_idx]
    
    cat(paste0("\n====================================================================\n"))
    cat(paste0("  Processing for Tau Distribution: Normal(μ=", current_threshold_mean, ", σ=", current_tau_sd, ")\n"))
    cat(paste0("   (Mean Iteration ", mean_idx, " of ", total_mean_iterations, " for current SD)\n"))
    cat(paste0("====================================================================\n"))
    
    cl <- makeCluster(NUM_CORES_TO_USE, type = "FORK")
    registerDoParallel(cl)
    cat(paste0("    Registered ", NUM_CORES_TO_USE, " parallel workers.\n"))
    cat(paste0("    Starting ", NUM_SEED_RUNS_TOTAL, " simulation runs (network+seed pairs)...\n"))
    
    # Progress tracking for the current (mean, sd) combination
    pb_env <- new.env()
    pb_env$completed_runs <- 0
    pb_env$total_runs <- NUM_SEED_RUNS_TOTAL
    
    list_of_results_for_this_mean_sd_combo <- foreach(
      run_idx = 1:NUM_SEED_RUNS_TOTAL,
      .combine = 'list',
      .multicombine = TRUE,
      .packages = c('igraph', 'dplyr', 'readr', 'intergraph', 'cluster'),
      .export = c('sweep_homoph_parameter', 'get_complex_plot',
                  'NETWORKS_DIR', 'current_threshold_mean', 'current_tau_sd', 
                  'IUL_VALUES_SWEEP', 'H_VALUES_SWEEP', 'SEEDING_STRATEGY_FIXED',
                  'NUM_NETWORK_INSTANCES_AVAILABLE', # For modulo network index
                  'N_NODES_GLOBAL' # Define if used, or ensure functions get N_nodes_arg
      ), 
      .errorhandling = 'pass',
      .final = function(x) { # .final to update progress after all workers return for this combo
        # This .final is for the whole foreach, not per worker iteration.
        # For per-worker progress, it's more complex with doParallel.
        # We will print progress after the parallel block.
        return(x)
      }
    ) %dopar% {
      
      # Determine network file index, cycling if NUM_SEED_RUNS_TOTAL > NUM_NETWORK_INSTANCES_AVAILABLE
      network_file_idx <- ((run_idx - 1) %% NUM_NETWORK_INSTANCES_AVAILABLE) + 1
      current_network_path <- paste0(NETWORKS_DIR, "ATP_network_simulated_1000_mur_", sprintf("%03d", network_file_idx), ".rds")
      
      # Minimal worker output to avoid clutter, can be enabled for debugging
      # worker_message_prefix <- paste0("      Worker ", Sys.getpid(), " [M:", sprintf("%.2f", current_threshold_mean), ",SD:", sprintf("%.2f", current_tau_sd), "] (Run ", run_idx, "): ")
      # cat(paste0(worker_message_prefix, "Net:", network_file_idx, "\n"))
      
      
      if (!file.exists(current_network_path)) return(NULL) 
      
      graph_for_this_run_ergm <- readRDS(current_network_path)
      graph_for_this_run <- asIgraph(graph_for_this_run_ergm)
      N_NODES_SPECIFIC_GRAPH <- vcount(graph_for_this_run)
      
      set.seed(run_idx * 3000 + round(current_threshold_mean * 100) + round(current_tau_sd * 100)) # Ajustar si current_tau_sd no está aquí
      node_mur_q_specific <- V(graph_for_this_run)$q_i
      # epsilon_q_i <- runif(N_NODES_SPECIFIC_GRAPH, min = -0.005, max = 0.005)
      # node_mur_q_specific <- node_mur_q_specific + epsilon_q_i
      # # Asegurar que q_i ruidoso esté entre 0 y 1
      # node_mur_q_specific[node_mur_q_specific < 0] <- 0
      # node_mur_q_specific[node_mur_q_specific > 1] <- 1
      
      node_degrees_specific <- igraph::degree(graph_for_this_run)
      
      attributes_for_distance_specific <- data.frame(
        age = V(graph_for_this_run)$age, educ_num = V(graph_for_this_run)$educ_num,
        race = as.factor(V(graph_for_this_run)$race), relig = as.factor(V(graph_for_this_run)$relig),
        sex = as.factor(V(graph_for_this_run)$sex)
      )
      d_ij_matrix <- as.matrix(daisy(attributes_for_distance_specific, metric = "gower"))
      
      # set.seed(run_idx * 4000 + round(current_threshold_mean * 100) + round(current_tau_sd * 100)) # Semilla diferente
      # epsilon_matrix_d_ij <- matrix(0, nrow = N_NODES_SPECIFIC_GRAPH, ncol = N_NODES_SPECIFIC_GRAPH)
      # upper_triangle_indices <- upper.tri(epsilon_matrix_d_ij)
      # noise_values_for_upper_triangle <- runif(sum(upper_triangle_indices), min = -0.025, max = 0.025)
      # epsilon_matrix_d_ij[upper_triangle_indices] <- noise_values_for_upper_triangle
      # epsilon_matrix_d_ij <- epsilon_matrix_d_ij + t(epsilon_matrix_d_ij) # Hacerla simétrica
      # 
      # d_ij_matrix <- d_ij_matrix + epsilon_matrix_d_ij
      # 
      # # Asegurar que d_ij ruidoso esté entre 0 y 1 (o el rango máximo de Gower, que es 1)
      # d_ij_matrix[d_ij_matrix < 0] <- 0
      # d_ij_matrix[d_ij_matrix > 1] <- 1
      # diag(d_ij_matrix) <- 0 # Asegurar que la diagonal sea cero
      
      # Thresholds
      set.seed(run_idx * 1000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000)) 
      node_thresholds_tau_frac_specific <- rnorm(
        n = N_NODES_SPECIFIC_GRAPH, mean = current_threshold_mean, sd = current_tau_sd
      )
      node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific < 0] <- 0
      node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific > 1] <- 1
      
      node_thresholds_count_for_cluster_specific <- round(node_thresholds_tau_frac_specific * node_degrees_specific)
      node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific > 0] <- 1
      node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific == 0] <- 0
      
      # if (SEEDING_STRATEGY_FIXED == "central") {
      #   degrees_current_graph <- igraph::degree(graph_for_this_run)
      #   primary_seed_for_this_run <- which.max(degrees_current_graph) 
      #   if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      # } else if (SEEDING_STRATEGY_FIXED == "random") { # Fallback a random si SEEDING_STRATEGY_FIXED no es "highest_degree"
      #   set.seed(run_idx * 2000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000))
      #   primary_seed_for_this_run <- as.integer(sample(V(graph_for_this_run), 1))
      # }
      # 
      if (SEEDING_STRATEGY_FIXED == "central") {
        primary_seed_for_this_run <- which.max(igraph::degree(graph_for_this_run))
        if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      } else if (SEEDING_STRATEGY_FIXED == "marginal") {
        lowest_10_percent <- sort(igraph::degree(graph_for_this_run), index.return = TRUE)$ix[1:ceiling(N_NODES_GLOBAL * 0.1)]
        set.seed(run_idx * 2000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000))
        primary_seed_for_this_run <- as.integer(sample(lowest_10_percent, 1))
      } else if (SEEDING_STRATEGY_FIXED == "closeness") {
        closeness_current_graph <- igraph::closeness(graph_for_this_run)
        primary_seed_for_this_run <- which.max(closeness_current_graph)
        if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      }  else if (SEEDING_STRATEGY_FIXED == "eigen") {
        eigen_scores <- igraph::eigen_centrality(graph_for_this_run)$vector
        primary_seed_for_this_run <- which.max(eigen_scores)
        if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      } else if (SEEDING_STRATEGY_FIXED == "random") {
        set.seed(run_idx * 2000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000))
        primary_seed_for_this_run <- as.integer(sample(V(graph_for_this_run), 1))
      }
      
      num_seeds_for_initial_cluster <- node_thresholds_count_for_cluster_specific[primary_seed_for_this_run]
      num_seeds_for_initial_cluster <- min(num_seeds_for_initial_cluster, N_NODES_SPECIFIC_GRAPH, node_degrees_specific[primary_seed_for_this_run] + 1)
      if (num_seeds_for_initial_cluster < 1) num_seeds_for_initial_cluster <- 1
      
      initial_infectors_for_this_sim_run <- c(primary_seed_for_this_run)
      if (num_seeds_for_initial_cluster > 1) {
        # ... (same logic for creating initial_infectors_for_this_sim_run as before) ...
        neighbors_of_primary <- as.numeric(neighbors(graph_for_this_run, primary_seed_for_this_run, mode="total"))
        num_additional_needed <- num_seeds_for_initial_cluster - 1
        if (length(neighbors_of_primary) >= num_additional_needed) {
          initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run, sample(neighbors_of_primary, num_additional_needed))
        } else {
          initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run, neighbors_of_primary)
          still_needed_more <- num_seeds_for_initial_cluster - length(initial_infectors_for_this_sim_run)
          if (still_needed_more > 0) {
            potential_others <- setdiff(1:N_NODES_SPECIFIC_GRAPH, initial_infectors_for_this_sim_run)
            if (length(potential_others) > 0) {
              initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run, sample(potential_others, min(length(potential_others), still_needed_more)))
            }
          }
        }
      }
      initial_infectors_for_this_sim_run <- unique(initial_infectors_for_this_sim_run)
      
      df_one_full_run <- sweep_homoph_parameter(
        primary_seed_id_arg = primary_seed_for_this_run, 
        N_nodes_arg = N_NODES_SPECIFIC_GRAPH, graph_obj_arg = graph_for_this_run,
        node_individual_thresholds_tau_arg = node_thresholds_tau_frac_specific,
        node_mur_q_arg = node_mur_q_specific,
        all_innovation_iul_Gamma_values = IUL_VALUES_SWEEP,
        all_social_distance_h_values = H_VALUES_SWEEP,
        initial_infectors_vector_arg = initial_infectors_for_this_sim_run,
        d_ij_matrix = d_ij_matrix
      )
      
      df_one_full_run$run_id <- run_idx 
      df_one_full_run$network_instance_file_idx <- network_file_idx
      df_one_full_run$threshold_mean_param <- current_threshold_mean 
      df_one_full_run$threshold_sd_param <- current_tau_sd
      df_one_full_run$N_nodes_actual <- N_NODES_SPECIFIC_GRAPH
      
      return(df_one_full_run)
    } # End foreach for NUM_SEED_RUNS_TOTAL
    
    stopCluster(cl)
    cat(paste0("    Parallel simulations for (Mean τ = ", current_threshold_mean, ", SD τ = ", current_tau_sd, ") finished.\n"))
    cat(paste0("    Collected results for ", length(list_of_results_for_this_mean_sd_combo), " runs.\n"))
    
    valid_results_indices <- !sapply(list_of_results_for_this_mean_sd_combo, function(x) inherits(x, "simpleError") || is.null(x) || nrow(x)==0)
    if(sum(valid_results_indices) > 0) {
      current_mean_sd_raw_df <- bind_rows(list_of_results_for_this_mean_sd_combo[valid_results_indices])
    } else {
      current_mean_sd_raw_df <- data.frame() 
    }
    
    if (nrow(current_mean_sd_raw_df) == 0) {
      cat(paste0("    WARNING: No valid results for (Mean τ = ", current_threshold_mean, ", SD τ = ", current_tau_sd, "). Skipping.\n"))
      current_sd_all_means_raw_results[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- NULL # Store NULL to know it was attempted
      current_sd_all_means_heatmap_data[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- NULL
      next 
    }
    current_sd_all_means_raw_results[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- current_mean_sd_raw_df
    
    # --- Phase Transition Identification for this (Mean, SD) combo ---
    phase_transition_data <- current_mean_sd_raw_df %>%
      arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
      group_by(run_id, social_distance_h) %>% 
      mutate(
        num_adopters_prop = num_adopters / N_nodes_actual,
        jump_in_adoption = num_adopters_prop - lag(num_adopters_prop, default = 0)
      ) %>%
      mutate(is_phase_transition_step = ifelse(jump_in_adoption >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)) %>%
      ungroup()
    
    first_transition_points <- phase_transition_data %>%
      filter(is_phase_transition_step == 1) %>%
      group_by(run_id, social_distance_h) %>%
      summarise(first_transition_IUL = min(innovation_iul_Gamma), .groups = 'drop')
    
    heatmap_data_list_temp <- list()
    for (iul_val in IUL_VALUES_SWEEP) {
      for (h_val in H_VALUES_SWEEP) {
        count_transitions <- first_transition_points %>%
          filter(social_distance_h == h_val, first_transition_IUL <= iul_val) %>% nrow()
        prop_transitions <- count_transitions / NUM_SEED_RUNS_TOTAL 
        heatmap_data_list_temp[[length(heatmap_data_list_temp) + 1]] <- data.frame(
          innovation_iul_Gamma = iul_val, social_distance_h = h_val,
          proportion_with_transition = prop_transitions,
          tau_mean_param = current_threshold_mean # Keep for potential internal faceting if needed
        )
      }
    }
    current_sd_all_means_heatmap_data[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- bind_rows(heatmap_data_list_temp)
    cat(paste0("    Phase transition analysis complete for (Mean τ = ", current_threshold_mean, ", SD τ = ", current_tau_sd, ").\n"))
    
  } # End inner loop for THRESHOLD_MEAN_SWEEP_LIST
  
  # --- Combine heatmap data for the CURRENT SD (from all its means) ---
  # And add the SD parameter for faceting the final plot
  if (length(current_sd_all_means_heatmap_data) > 0) {
    valid_heatmap_data_for_sd <- current_sd_all_means_heatmap_data[!sapply(current_sd_all_means_heatmap_data, is.null)]
    if (length(valid_heatmap_data_for_sd) > 0) {
      df_heatmap_this_sd <- bind_rows(valid_heatmap_data_for_sd)
      df_heatmap_this_sd$tau_sd_param <- current_tau_sd 
      all_sds_heatmap_data_list[[paste0("sd_", sprintf("%.2f", current_tau_sd))]] <- df_heatmap_this_sd
      
      # --- Generate and Save Plot for CURRENT SD ---
      df_heatmap_this_sd$social_distance_h_factor <- factor(sprintf("%.2f", df_heatmap_this_sd$social_distance_h))
      df_heatmap_this_sd$tau_mean_facet_label <- factor(paste0("Mean τ = ", sprintf("%.2f", df_heatmap_this_sd$tau_mean_param)))
      
      plot_for_current_sd <- ggplot(df_heatmap_this_sd, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = proportion_with_transition)) +
        geom_tile(color = "white", lwd = 0.15) +
        scale_fill_viridis_c(name = "Prop. of Runs\nwith Phase Transition\n(Jump >= 0.5)", limits = c(0, 1), option="plasma", n.breaks=6) +
        facet_wrap(~tau_mean_facet_label, ncol = 2) + 
        labs(
          x = expression(paste("Intrinsic Innovation Utility (", Gamma, ")")),
          y = "Maximum Social Distance (h)",
          title = paste("Phase Transition Maps for", CURRENT_GRAPH_TYPE_LABEL, "- Tau SD =", sprintf("%.2f", current_tau_sd)),
          subtitle = paste("Thresholds ~ N(μ=variable, σ=", sprintf("%.2f", current_tau_sd), "), ", NUM_SEED_RUNS_TOTAL, " runs per (Γ,h) per panel", sep="")
        ) +
        scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0,0)) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size=12),
          plot.subtitle = element_text(hjust = 0.5, size=8),
          axis.text.x = element_text(angle = 45, hjust = 1, size=7),
          axis.text.y = element_text(size=7), legend.position = "right",
          strip.text = element_text(face="bold", size=9) 
        )
      
      plot_filename_current_sd <- paste0(PLOTS_DIR, "phase_transition_sd", sprintf("%.2f", current_tau_sd), "_means03-06_", SEEDING_STRATEGY_FIXED, ".pdf")
      ggsave(plot_filename_current_sd, plot_for_current_sd, width = 10, height = 8) # Adjusted for potentially 4 means
      cat(paste0("    Saved plot: ", plot_filename_current_sd, "\n"))
    }
  }
  
  # --- Save Raw Data for CURRENT SD ---
  raw_data_filename_current_sd <- paste0(RESULTS_DIR, "phase_transition_raw_sd", sprintf("%.2f", current_tau_sd), "_means03-06_", SEEDING_STRATEGY_FIXED, ".rds")
  saveRDS(current_sd_all_means_raw_results, raw_data_filename_current_sd) # Save list of DFs (one per mean)
  cat(paste0("    Saved raw data: ", raw_data_filename_current_sd, "\n"))
  
  # Store the raw results for this SD in the grand list
  all_sds_raw_results_list[[paste0("sd_", sprintf("%.2f", current_tau_sd))]] <- current_sd_all_means_raw_results
  
} # End outer loop for TAU_NORMAL_SD_SWEEP_LIST
time_fin <- Sys.time()
time_total_parallel_1 <- difftime(time_fin, time_init, units = "auto")

cat("\nAll simulation sweeps complete.\n")

# -----------------------------------------------------------------------------
# 6. Final Grand Data Saving (optional, as data is saved per SD)
# -----------------------------------------------------------------------------
cat("Saving grand combined data objects...\n")

# Grand combined heatmap data (all SDs and all Means)
grand_combined_heatmap_df <- bind_rows(all_sds_heatmap_data_list)
saveRDS(grand_combined_heatmap_df, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_heatmap_data_all_sds_means_", SEEDING_STRATEGY_FIXED, ".rds"))

# Grand combined raw results (list of lists)
saveRDS(all_sds_raw_results_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_raw_results_all_sds_means_", SEEDING_STRATEGY_FIXED, ".rds"))

cat("All data saved. Script finished.\n")









# -----------------------------------------------------------------------------
# 3. Main Simulation Loop (Outer loop for SD, Inner loop for Mean)
# -----------------------------------------------------------------------------
cat("Starting grand simulation sweep...\n")

SEEDING_STRATEGY_FIXED <- "closeness" #"marginal" #"eigen" #"central" #"random"

# This will store all heatmap dataframes, one for each SD
all_sds_heatmap_data_list <- list()
# This will store all raw results, one list element per SD, 
# where each element is itself a list of results per Mean
all_sds_raw_results_list <- list()

total_sd_iterations <- length(TAU_NORMAL_SD_SWEEP_LIST)

time_init <- Sys.time()
for (sd_idx in 1:total_sd_iterations) {
  current_tau_sd <- TAU_NORMAL_SD_SWEEP_LIST[sd_idx]
  
  cat(paste0("\n####################################################################\n"))
  cat(paste0("Processing for Tau Distribution Standard Deviation (σ): ", current_tau_sd, "\n"))
  cat(paste0(" (SD Iteration ", sd_idx, " of ", total_sd_iterations, ")\n"))
  cat(paste0("####################################################################\n"))
  
  # This will store raw results for the current SD, one list element per Mean
  current_sd_all_means_raw_results <- list()
  # This will store heatmap data for the current SD, one list element per Mean (for faceting)
  current_sd_all_means_heatmap_data <- list()
  
  total_mean_iterations <- length(THRESHOLD_MEAN_SWEEP_LIST)
  for (mean_idx in 1:total_mean_iterations) {
    current_threshold_mean <- THRESHOLD_MEAN_SWEEP_LIST[mean_idx]
    
    cat(paste0("\n====================================================================\n"))
    cat(paste0("  Processing for Tau Distribution: Normal(μ=", current_threshold_mean, ", σ=", current_tau_sd, ")\n"))
    cat(paste0("   (Mean Iteration ", mean_idx, " of ", total_mean_iterations, " for current SD)\n"))
    cat(paste0("====================================================================\n"))
    
    cl <- makeCluster(NUM_CORES_TO_USE, type = "FORK")
    registerDoParallel(cl)
    cat(paste0("    Registered ", NUM_CORES_TO_USE, " parallel workers.\n"))
    cat(paste0("    Starting ", NUM_SEED_RUNS_TOTAL, " simulation runs (network+seed pairs)...\n"))
    
    # Progress tracking for the current (mean, sd) combination
    pb_env <- new.env()
    pb_env$completed_runs <- 0
    pb_env$total_runs <- NUM_SEED_RUNS_TOTAL
    
    list_of_results_for_this_mean_sd_combo <- foreach(
      run_idx = 1:NUM_SEED_RUNS_TOTAL,
      .combine = 'list',
      .multicombine = TRUE,
      .packages = c('igraph', 'dplyr', 'readr', 'intergraph', 'cluster'),
      .export = c('sweep_homoph_parameter', 'get_complex_plot',
                  'NETWORKS_DIR', 'current_threshold_mean', 'current_tau_sd', 
                  'IUL_VALUES_SWEEP', 'H_VALUES_SWEEP', 'SEEDING_STRATEGY_FIXED',
                  'NUM_NETWORK_INSTANCES_AVAILABLE', # For modulo network index
                  'N_NODES_GLOBAL' # Define if used, or ensure functions get N_nodes_arg
      ), 
      .errorhandling = 'pass',
      .final = function(x) { # .final to update progress after all workers return for this combo
        # This .final is for the whole foreach, not per worker iteration.
        # For per-worker progress, it's more complex with doParallel.
        # We will print progress after the parallel block.
        return(x)
      }
    ) %dopar% {
      
      # Determine network file index, cycling if NUM_SEED_RUNS_TOTAL > NUM_NETWORK_INSTANCES_AVAILABLE
      network_file_idx <- ((run_idx - 1) %% NUM_NETWORK_INSTANCES_AVAILABLE) + 1
      current_network_path <- paste0(NETWORKS_DIR, "ATP_network_simulated_1000_mur_", sprintf("%03d", network_file_idx), ".rds")
      
      # Minimal worker output to avoid clutter, can be enabled for debugging
      # worker_message_prefix <- paste0("      Worker ", Sys.getpid(), " [M:", sprintf("%.2f", current_threshold_mean), ",SD:", sprintf("%.2f", current_tau_sd), "] (Run ", run_idx, "): ")
      # cat(paste0(worker_message_prefix, "Net:", network_file_idx, "\n"))
      
      
      if (!file.exists(current_network_path)) return(NULL) 
      
      graph_for_this_run_ergm <- readRDS(current_network_path)
      graph_for_this_run <- asIgraph(graph_for_this_run_ergm)
      N_NODES_SPECIFIC_GRAPH <- vcount(graph_for_this_run)
      
      set.seed(run_idx * 3000 + round(current_threshold_mean * 100) + round(current_tau_sd * 100)) # Ajustar si current_tau_sd no está aquí
      node_mur_q_specific <- V(graph_for_this_run)$q_i
      # epsilon_q_i <- runif(N_NODES_SPECIFIC_GRAPH, min = -0.005, max = 0.005)
      # node_mur_q_specific <- node_mur_q_specific + epsilon_q_i
      # # Asegurar que q_i ruidoso esté entre 0 y 1
      # node_mur_q_specific[node_mur_q_specific < 0] <- 0
      # node_mur_q_specific[node_mur_q_specific > 1] <- 1
      
      node_degrees_specific <- igraph::degree(graph_for_this_run)
      
      attributes_for_distance_specific <- data.frame(
        age = V(graph_for_this_run)$age, educ_num = V(graph_for_this_run)$educ_num,
        race = as.factor(V(graph_for_this_run)$race), relig = as.factor(V(graph_for_this_run)$relig),
        sex = as.factor(V(graph_for_this_run)$sex)
      )
      d_ij_matrix <- as.matrix(daisy(attributes_for_distance_specific, metric = "gower"))
      
      # set.seed(run_idx * 4000 + round(current_threshold_mean * 100) + round(current_tau_sd * 100)) # Semilla diferente
      # epsilon_matrix_d_ij <- matrix(0, nrow = N_NODES_SPECIFIC_GRAPH, ncol = N_NODES_SPECIFIC_GRAPH)
      # upper_triangle_indices <- upper.tri(epsilon_matrix_d_ij)
      # noise_values_for_upper_triangle <- runif(sum(upper_triangle_indices), min = -0.025, max = 0.025)
      # epsilon_matrix_d_ij[upper_triangle_indices] <- noise_values_for_upper_triangle
      # epsilon_matrix_d_ij <- epsilon_matrix_d_ij + t(epsilon_matrix_d_ij) # Hacerla simétrica
      # 
      # d_ij_matrix <- d_ij_matrix + epsilon_matrix_d_ij
      # 
      # # Asegurar que d_ij ruidoso esté entre 0 y 1 (o el rango máximo de Gower, que es 1)
      # d_ij_matrix[d_ij_matrix < 0] <- 0
      # d_ij_matrix[d_ij_matrix > 1] <- 1
      # diag(d_ij_matrix) <- 0 # Asegurar que la diagonal sea cero
      
      # Thresholds
      set.seed(run_idx * 1000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000)) 
      node_thresholds_tau_frac_specific <- rnorm(
        n = N_NODES_SPECIFIC_GRAPH, mean = current_threshold_mean, sd = current_tau_sd
      )
      node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific < 0] <- 0
      node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific > 1] <- 1
      
      node_thresholds_count_for_cluster_specific <- round(node_thresholds_tau_frac_specific * node_degrees_specific)
      node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific > 0] <- 1
      node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific == 0] <- 0
      
      # if (SEEDING_STRATEGY_FIXED == "central") {
      #   degrees_current_graph <- igraph::degree(graph_for_this_run)
      #   primary_seed_for_this_run <- which.max(degrees_current_graph) 
      #   if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      # } else if (SEEDING_STRATEGY_FIXED == "random") { # Fallback a random si SEEDING_STRATEGY_FIXED no es "highest_degree"
      #   set.seed(run_idx * 2000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000))
      #   primary_seed_for_this_run <- as.integer(sample(V(graph_for_this_run), 1))
      # }
      # 
      if (SEEDING_STRATEGY_FIXED == "central") {
        primary_seed_for_this_run <- which.max(igraph::degree(graph_for_this_run))
        if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      } else if (SEEDING_STRATEGY_FIXED == "marginal") {
        lowest_10_percent <- sort(igraph::degree(graph_for_this_run), index.return = TRUE)$ix[1:ceiling(N_NODES_GLOBAL * 0.1)]
        set.seed(run_idx * 2000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000))
        primary_seed_for_this_run <- as.integer(sample(lowest_10_percent, 1))
      } else if (SEEDING_STRATEGY_FIXED == "closeness") {
        closeness_current_graph <- igraph::closeness(graph_for_this_run)
        primary_seed_for_this_run <- which.max(closeness_current_graph)
        if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      }  else if (SEEDING_STRATEGY_FIXED == "eigen") {
        eigen_scores <- igraph::eigen_centrality(graph_for_this_run)$vector
        primary_seed_for_this_run <- which.max(eigen_scores)
        if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      } else if (SEEDING_STRATEGY_FIXED == "random") {
        set.seed(run_idx * 2000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000))
        primary_seed_for_this_run <- as.integer(sample(V(graph_for_this_run), 1))
      }
      
      num_seeds_for_initial_cluster <- node_thresholds_count_for_cluster_specific[primary_seed_for_this_run]
      num_seeds_for_initial_cluster <- min(num_seeds_for_initial_cluster, N_NODES_SPECIFIC_GRAPH, node_degrees_specific[primary_seed_for_this_run] + 1)
      if (num_seeds_for_initial_cluster < 1) num_seeds_for_initial_cluster <- 1
      
      initial_infectors_for_this_sim_run <- c(primary_seed_for_this_run)
      if (num_seeds_for_initial_cluster > 1) {
        # ... (same logic for creating initial_infectors_for_this_sim_run as before) ...
        neighbors_of_primary <- as.numeric(neighbors(graph_for_this_run, primary_seed_for_this_run, mode="total"))
        num_additional_needed <- num_seeds_for_initial_cluster - 1
        if (length(neighbors_of_primary) >= num_additional_needed) {
          initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run, sample(neighbors_of_primary, num_additional_needed))
        } else {
          initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run, neighbors_of_primary)
          still_needed_more <- num_seeds_for_initial_cluster - length(initial_infectors_for_this_sim_run)
          if (still_needed_more > 0) {
            potential_others <- setdiff(1:N_NODES_SPECIFIC_GRAPH, initial_infectors_for_this_sim_run)
            if (length(potential_others) > 0) {
              initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run, sample(potential_others, min(length(potential_others), still_needed_more)))
            }
          }
        }
      }
      initial_infectors_for_this_sim_run <- unique(initial_infectors_for_this_sim_run)
      
      df_one_full_run <- sweep_homoph_parameter(
        primary_seed_id_arg = primary_seed_for_this_run, 
        N_nodes_arg = N_NODES_SPECIFIC_GRAPH, graph_obj_arg = graph_for_this_run,
        node_individual_thresholds_tau_arg = node_thresholds_tau_frac_specific,
        node_mur_q_arg = node_mur_q_specific,
        all_innovation_iul_Gamma_values = IUL_VALUES_SWEEP,
        all_social_distance_h_values = H_VALUES_SWEEP,
        initial_infectors_vector_arg = initial_infectors_for_this_sim_run,
        d_ij_matrix = d_ij_matrix
      )
      
      df_one_full_run$run_id <- run_idx 
      df_one_full_run$network_instance_file_idx <- network_file_idx
      df_one_full_run$threshold_mean_param <- current_threshold_mean 
      df_one_full_run$threshold_sd_param <- current_tau_sd
      df_one_full_run$N_nodes_actual <- N_NODES_SPECIFIC_GRAPH
      
      return(df_one_full_run)
    } # End foreach for NUM_SEED_RUNS_TOTAL
    
    stopCluster(cl)
    cat(paste0("    Parallel simulations for (Mean τ = ", current_threshold_mean, ", SD τ = ", current_tau_sd, ") finished.\n"))
    cat(paste0("    Collected results for ", length(list_of_results_for_this_mean_sd_combo), " runs.\n"))
    
    valid_results_indices <- !sapply(list_of_results_for_this_mean_sd_combo, function(x) inherits(x, "simpleError") || is.null(x) || nrow(x)==0)
    if(sum(valid_results_indices) > 0) {
      current_mean_sd_raw_df <- bind_rows(list_of_results_for_this_mean_sd_combo[valid_results_indices])
    } else {
      current_mean_sd_raw_df <- data.frame() 
    }
    
    if (nrow(current_mean_sd_raw_df) == 0) {
      cat(paste0("    WARNING: No valid results for (Mean τ = ", current_threshold_mean, ", SD τ = ", current_tau_sd, "). Skipping.\n"))
      current_sd_all_means_raw_results[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- NULL # Store NULL to know it was attempted
      current_sd_all_means_heatmap_data[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- NULL
      next 
    }
    current_sd_all_means_raw_results[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- current_mean_sd_raw_df
    
    # --- Phase Transition Identification for this (Mean, SD) combo ---
    phase_transition_data <- current_mean_sd_raw_df %>%
      arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
      group_by(run_id, social_distance_h) %>% 
      mutate(
        num_adopters_prop = num_adopters / N_nodes_actual,
        jump_in_adoption = num_adopters_prop - lag(num_adopters_prop, default = 0)
      ) %>%
      mutate(is_phase_transition_step = ifelse(jump_in_adoption >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)) %>%
      ungroup()
    
    first_transition_points <- phase_transition_data %>%
      filter(is_phase_transition_step == 1) %>%
      group_by(run_id, social_distance_h) %>%
      summarise(first_transition_IUL = min(innovation_iul_Gamma), .groups = 'drop')
    
    heatmap_data_list_temp <- list()
    for (iul_val in IUL_VALUES_SWEEP) {
      for (h_val in H_VALUES_SWEEP) {
        count_transitions <- first_transition_points %>%
          filter(social_distance_h == h_val, first_transition_IUL <= iul_val) %>% nrow()
        prop_transitions <- count_transitions / NUM_SEED_RUNS_TOTAL 
        heatmap_data_list_temp[[length(heatmap_data_list_temp) + 1]] <- data.frame(
          innovation_iul_Gamma = iul_val, social_distance_h = h_val,
          proportion_with_transition = prop_transitions,
          tau_mean_param = current_threshold_mean # Keep for potential internal faceting if needed
        )
      }
    }
    current_sd_all_means_heatmap_data[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- bind_rows(heatmap_data_list_temp)
    cat(paste0("    Phase transition analysis complete for (Mean τ = ", current_threshold_mean, ", SD τ = ", current_tau_sd, ").\n"))
    
  } # End inner loop for THRESHOLD_MEAN_SWEEP_LIST
  
  # --- Combine heatmap data for the CURRENT SD (from all its means) ---
  # And add the SD parameter for faceting the final plot
  if (length(current_sd_all_means_heatmap_data) > 0) {
    valid_heatmap_data_for_sd <- current_sd_all_means_heatmap_data[!sapply(current_sd_all_means_heatmap_data, is.null)]
    if (length(valid_heatmap_data_for_sd) > 0) {
      df_heatmap_this_sd <- bind_rows(valid_heatmap_data_for_sd)
      df_heatmap_this_sd$tau_sd_param <- current_tau_sd 
      all_sds_heatmap_data_list[[paste0("sd_", sprintf("%.2f", current_tau_sd))]] <- df_heatmap_this_sd
      
      # --- Generate and Save Plot for CURRENT SD ---
      df_heatmap_this_sd$social_distance_h_factor <- factor(sprintf("%.2f", df_heatmap_this_sd$social_distance_h))
      df_heatmap_this_sd$tau_mean_facet_label <- factor(paste0("Mean τ = ", sprintf("%.2f", df_heatmap_this_sd$tau_mean_param)))
      
      plot_for_current_sd <- ggplot(df_heatmap_this_sd, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = proportion_with_transition)) +
        geom_tile(color = "white", lwd = 0.15) +
        scale_fill_viridis_c(name = "Prop. of Runs\nwith Phase Transition\n(Jump >= 0.5)", limits = c(0, 1), option="plasma", n.breaks=6) +
        facet_wrap(~tau_mean_facet_label, ncol = 2) + 
        labs(
          x = expression(paste("Intrinsic Innovation Utility (", Gamma, ")")),
          y = "Maximum Social Distance (h)",
          title = paste("Phase Transition Maps for", CURRENT_GRAPH_TYPE_LABEL, "- Tau SD =", sprintf("%.2f", current_tau_sd)),
          subtitle = paste("Thresholds ~ N(μ=variable, σ=", sprintf("%.2f", current_tau_sd), "), ", NUM_SEED_RUNS_TOTAL, " runs per (Γ,h) per panel", sep="")
        ) +
        scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0,0)) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size=12),
          plot.subtitle = element_text(hjust = 0.5, size=8),
          axis.text.x = element_text(angle = 45, hjust = 1, size=7),
          axis.text.y = element_text(size=7), legend.position = "right",
          strip.text = element_text(face="bold", size=9) 
        )
      
      plot_filename_current_sd <- paste0(PLOTS_DIR, "phase_transition_sd", sprintf("%.2f", current_tau_sd), "_means03-06_", SEEDING_STRATEGY_FIXED, ".pdf")
      ggsave(plot_filename_current_sd, plot_for_current_sd, width = 10, height = 8) # Adjusted for potentially 4 means
      cat(paste0("    Saved plot: ", plot_filename_current_sd, "\n"))
    }
  }
  
  # --- Save Raw Data for CURRENT SD ---
  raw_data_filename_current_sd <- paste0(RESULTS_DIR, "phase_transition_raw_sd", sprintf("%.2f", current_tau_sd), "_means03-06_", SEEDING_STRATEGY_FIXED, ".rds")
  saveRDS(current_sd_all_means_raw_results, raw_data_filename_current_sd) # Save list of DFs (one per mean)
  cat(paste0("    Saved raw data: ", raw_data_filename_current_sd, "\n"))
  
  # Store the raw results for this SD in the grand list
  all_sds_raw_results_list[[paste0("sd_", sprintf("%.2f", current_tau_sd))]] <- current_sd_all_means_raw_results
  
} # End outer loop for TAU_NORMAL_SD_SWEEP_LIST
time_fin <- Sys.time()
time_total_parallel_2 <- difftime(time_fin, time_init, units = "auto")

cat("\nAll simulation sweeps complete.\n")

# -----------------------------------------------------------------------------
# 6. Final Grand Data Saving (optional, as data is saved per SD)
# -----------------------------------------------------------------------------
cat("Saving grand combined data objects...\n")

# Grand combined heatmap data (all SDs and all Means)
grand_combined_heatmap_df <- bind_rows(all_sds_heatmap_data_list)
saveRDS(grand_combined_heatmap_df, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_heatmap_data_all_sds_means_", SEEDING_STRATEGY_FIXED, ".rds"))

# Grand combined raw results (list of lists)
saveRDS(all_sds_raw_results_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_raw_results_all_sds_means_", SEEDING_STRATEGY_FIXED, ".rds"))

cat("All data saved. Script finished.\n")














# -----------------------------------------------------------------------------
# 3. Main Simulation Loop (Outer loop for SD, Inner loop for Mean)
# -----------------------------------------------------------------------------
cat("Starting grand simulation sweep...\n")

SEEDING_STRATEGY_FIXED <- "marginal" #"eigen" #"central" #"random"

# This will store all heatmap dataframes, one for each SD
all_sds_heatmap_data_list <- list()
# This will store all raw results, one list element per SD, 
# where each element is itself a list of results per Mean
all_sds_raw_results_list <- list()

total_sd_iterations <- length(TAU_NORMAL_SD_SWEEP_LIST)

time_init <- Sys.time()
for (sd_idx in 1:total_sd_iterations) {
  current_tau_sd <- TAU_NORMAL_SD_SWEEP_LIST[sd_idx]
  
  cat(paste0("\n####################################################################\n"))
  cat(paste0("Processing for Tau Distribution Standard Deviation (σ): ", current_tau_sd, "\n"))
  cat(paste0(" (SD Iteration ", sd_idx, " of ", total_sd_iterations, ")\n"))
  cat(paste0("####################################################################\n"))
  
  # This will store raw results for the current SD, one list element per Mean
  current_sd_all_means_raw_results <- list()
  # This will store heatmap data for the current SD, one list element per Mean (for faceting)
  current_sd_all_means_heatmap_data <- list()
  
  total_mean_iterations <- length(THRESHOLD_MEAN_SWEEP_LIST)
  for (mean_idx in 1:total_mean_iterations) {
    current_threshold_mean <- THRESHOLD_MEAN_SWEEP_LIST[mean_idx]
    
    cat(paste0("\n====================================================================\n"))
    cat(paste0("  Processing for Tau Distribution: Normal(μ=", current_threshold_mean, ", σ=", current_tau_sd, ")\n"))
    cat(paste0("   (Mean Iteration ", mean_idx, " of ", total_mean_iterations, " for current SD)\n"))
    cat(paste0("====================================================================\n"))
    
    cl <- makeCluster(NUM_CORES_TO_USE, type = "FORK")
    registerDoParallel(cl)
    cat(paste0("    Registered ", NUM_CORES_TO_USE, " parallel workers.\n"))
    cat(paste0("    Starting ", NUM_SEED_RUNS_TOTAL, " simulation runs (network+seed pairs)...\n"))
    
    # Progress tracking for the current (mean, sd) combination
    pb_env <- new.env()
    pb_env$completed_runs <- 0
    pb_env$total_runs <- NUM_SEED_RUNS_TOTAL
    
    list_of_results_for_this_mean_sd_combo <- foreach(
      run_idx = 1:NUM_SEED_RUNS_TOTAL,
      .combine = 'list',
      .multicombine = TRUE,
      .packages = c('igraph', 'dplyr', 'readr', 'intergraph', 'cluster'),
      .export = c('sweep_homoph_parameter', 'get_complex_plot',
                  'NETWORKS_DIR', 'current_threshold_mean', 'current_tau_sd', 
                  'IUL_VALUES_SWEEP', 'H_VALUES_SWEEP', 'SEEDING_STRATEGY_FIXED',
                  'NUM_NETWORK_INSTANCES_AVAILABLE', # For modulo network index
                  'N_NODES_GLOBAL' # Define if used, or ensure functions get N_nodes_arg
      ), 
      .errorhandling = 'pass',
      .final = function(x) { # .final to update progress after all workers return for this combo
        # This .final is for the whole foreach, not per worker iteration.
        # For per-worker progress, it's more complex with doParallel.
        # We will print progress after the parallel block.
        return(x)
      }
    ) %dopar% {
      
      # Determine network file index, cycling if NUM_SEED_RUNS_TOTAL > NUM_NETWORK_INSTANCES_AVAILABLE
      network_file_idx <- ((run_idx - 1) %% NUM_NETWORK_INSTANCES_AVAILABLE) + 1
      current_network_path <- paste0(NETWORKS_DIR, "ATP_network_simulated_1000_mur_", sprintf("%03d", network_file_idx), ".rds")
      
      # Minimal worker output to avoid clutter, can be enabled for debugging
      # worker_message_prefix <- paste0("      Worker ", Sys.getpid(), " [M:", sprintf("%.2f", current_threshold_mean), ",SD:", sprintf("%.2f", current_tau_sd), "] (Run ", run_idx, "): ")
      # cat(paste0(worker_message_prefix, "Net:", network_file_idx, "\n"))
      
      
      if (!file.exists(current_network_path)) return(NULL) 
      
      graph_for_this_run_ergm <- readRDS(current_network_path)
      graph_for_this_run <- asIgraph(graph_for_this_run_ergm)
      N_NODES_SPECIFIC_GRAPH <- vcount(graph_for_this_run)
      
      set.seed(run_idx * 3000 + round(current_threshold_mean * 100) + round(current_tau_sd * 100)) # Ajustar si current_tau_sd no está aquí
      node_mur_q_specific <- V(graph_for_this_run)$q_i
      # epsilon_q_i <- runif(N_NODES_SPECIFIC_GRAPH, min = -0.005, max = 0.005)
      # node_mur_q_specific <- node_mur_q_specific + epsilon_q_i
      # # Asegurar que q_i ruidoso esté entre 0 y 1
      # node_mur_q_specific[node_mur_q_specific < 0] <- 0
      # node_mur_q_specific[node_mur_q_specific > 1] <- 1
      
      node_degrees_specific <- igraph::degree(graph_for_this_run)
      
      attributes_for_distance_specific <- data.frame(
        age = V(graph_for_this_run)$age, educ_num = V(graph_for_this_run)$educ_num,
        race = as.factor(V(graph_for_this_run)$race), relig = as.factor(V(graph_for_this_run)$relig),
        sex = as.factor(V(graph_for_this_run)$sex)
      )
      d_ij_matrix <- as.matrix(daisy(attributes_for_distance_specific, metric = "gower"))
      
      # set.seed(run_idx * 4000 + round(current_threshold_mean * 100) + round(current_tau_sd * 100)) # Semilla diferente
      # epsilon_matrix_d_ij <- matrix(0, nrow = N_NODES_SPECIFIC_GRAPH, ncol = N_NODES_SPECIFIC_GRAPH)
      # upper_triangle_indices <- upper.tri(epsilon_matrix_d_ij)
      # noise_values_for_upper_triangle <- runif(sum(upper_triangle_indices), min = -0.025, max = 0.025)
      # epsilon_matrix_d_ij[upper_triangle_indices] <- noise_values_for_upper_triangle
      # epsilon_matrix_d_ij <- epsilon_matrix_d_ij + t(epsilon_matrix_d_ij) # Hacerla simétrica
      # 
      # d_ij_matrix <- d_ij_matrix + epsilon_matrix_d_ij
      # 
      # # Asegurar que d_ij ruidoso esté entre 0 y 1 (o el rango máximo de Gower, que es 1)
      # d_ij_matrix[d_ij_matrix < 0] <- 0
      # d_ij_matrix[d_ij_matrix > 1] <- 1
      # diag(d_ij_matrix) <- 0 # Asegurar que la diagonal sea cero
      
      # Thresholds
      set.seed(run_idx * 1000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000)) 
      node_thresholds_tau_frac_specific <- rnorm(
        n = N_NODES_SPECIFIC_GRAPH, mean = current_threshold_mean, sd = current_tau_sd
      )
      node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific < 0] <- 0
      node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific > 1] <- 1
      
      node_thresholds_count_for_cluster_specific <- round(node_thresholds_tau_frac_specific * node_degrees_specific)
      node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific > 0] <- 1
      node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific == 0] <- 0
      
      # if (SEEDING_STRATEGY_FIXED == "central") {
      #   degrees_current_graph <- igraph::degree(graph_for_this_run)
      #   primary_seed_for_this_run <- which.max(degrees_current_graph) 
      #   if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      # } else if (SEEDING_STRATEGY_FIXED == "random") { # Fallback a random si SEEDING_STRATEGY_FIXED no es "highest_degree"
      #   set.seed(run_idx * 2000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000))
      #   primary_seed_for_this_run <- as.integer(sample(V(graph_for_this_run), 1))
      # }
      # 
      if (SEEDING_STRATEGY_FIXED == "central") {
        primary_seed_for_this_run <- which.max(igraph::degree(graph_for_this_run))
        if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      } else if (SEEDING_STRATEGY_FIXED == "marginal") {
        lowest_10_percent <- sort(igraph::degree(graph_for_this_run), index.return = TRUE)$ix[1:ceiling(N_NODES_GLOBAL * 0.1)]
        set.seed(run_idx * 2000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000))
        primary_seed_for_this_run <- as.integer(sample(lowest_10_percent, 1))
      } else if (SEEDING_STRATEGY_FIXED == "closeness") {
        closeness_current_graph <- igraph::closeness(graph_for_this_run)
        primary_seed_for_this_run <- which.max(closeness_current_graph)
        if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      }  else if (SEEDING_STRATEGY_FIXED == "eigen") {
        eigen_scores <- igraph::eigen_centrality(graph_for_this_run)$vector
        primary_seed_for_this_run <- which.max(eigen_scores)
        if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      } else if (SEEDING_STRATEGY_FIXED == "random") {
        set.seed(run_idx * 2000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000))
        primary_seed_for_this_run <- as.integer(sample(V(graph_for_this_run), 1))
      }
      
      num_seeds_for_initial_cluster <- node_thresholds_count_for_cluster_specific[primary_seed_for_this_run]
      num_seeds_for_initial_cluster <- min(num_seeds_for_initial_cluster, N_NODES_SPECIFIC_GRAPH, node_degrees_specific[primary_seed_for_this_run] + 1)
      if (num_seeds_for_initial_cluster < 1) num_seeds_for_initial_cluster <- 1
      
      initial_infectors_for_this_sim_run <- c(primary_seed_for_this_run)
      if (num_seeds_for_initial_cluster > 1) {
        # ... (same logic for creating initial_infectors_for_this_sim_run as before) ...
        neighbors_of_primary <- as.numeric(neighbors(graph_for_this_run, primary_seed_for_this_run, mode="total"))
        num_additional_needed <- num_seeds_for_initial_cluster - 1
        if (length(neighbors_of_primary) >= num_additional_needed) {
          initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run, sample(neighbors_of_primary, num_additional_needed))
        } else {
          initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run, neighbors_of_primary)
          still_needed_more <- num_seeds_for_initial_cluster - length(initial_infectors_for_this_sim_run)
          if (still_needed_more > 0) {
            potential_others <- setdiff(1:N_NODES_SPECIFIC_GRAPH, initial_infectors_for_this_sim_run)
            if (length(potential_others) > 0) {
              initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run, sample(potential_others, min(length(potential_others), still_needed_more)))
            }
          }
        }
      }
      initial_infectors_for_this_sim_run <- unique(initial_infectors_for_this_sim_run)
      
      df_one_full_run <- sweep_homoph_parameter(
        primary_seed_id_arg = primary_seed_for_this_run, 
        N_nodes_arg = N_NODES_SPECIFIC_GRAPH, graph_obj_arg = graph_for_this_run,
        node_individual_thresholds_tau_arg = node_thresholds_tau_frac_specific,
        node_mur_q_arg = node_mur_q_specific,
        all_innovation_iul_Gamma_values = IUL_VALUES_SWEEP,
        all_social_distance_h_values = H_VALUES_SWEEP,
        initial_infectors_vector_arg = initial_infectors_for_this_sim_run,
        d_ij_matrix = d_ij_matrix
      )
      
      df_one_full_run$run_id <- run_idx 
      df_one_full_run$network_instance_file_idx <- network_file_idx
      df_one_full_run$threshold_mean_param <- current_threshold_mean 
      df_one_full_run$threshold_sd_param <- current_tau_sd
      df_one_full_run$N_nodes_actual <- N_NODES_SPECIFIC_GRAPH
      
      return(df_one_full_run)
    } # End foreach for NUM_SEED_RUNS_TOTAL
    
    stopCluster(cl)
    cat(paste0("    Parallel simulations for (Mean τ = ", current_threshold_mean, ", SD τ = ", current_tau_sd, ") finished.\n"))
    cat(paste0("    Collected results for ", length(list_of_results_for_this_mean_sd_combo), " runs.\n"))
    
    valid_results_indices <- !sapply(list_of_results_for_this_mean_sd_combo, function(x) inherits(x, "simpleError") || is.null(x) || nrow(x)==0)
    if(sum(valid_results_indices) > 0) {
      current_mean_sd_raw_df <- bind_rows(list_of_results_for_this_mean_sd_combo[valid_results_indices])
    } else {
      current_mean_sd_raw_df <- data.frame() 
    }
    
    if (nrow(current_mean_sd_raw_df) == 0) {
      cat(paste0("    WARNING: No valid results for (Mean τ = ", current_threshold_mean, ", SD τ = ", current_tau_sd, "). Skipping.\n"))
      current_sd_all_means_raw_results[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- NULL # Store NULL to know it was attempted
      current_sd_all_means_heatmap_data[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- NULL
      next 
    }
    current_sd_all_means_raw_results[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- current_mean_sd_raw_df
    
    # --- Phase Transition Identification for this (Mean, SD) combo ---
    phase_transition_data <- current_mean_sd_raw_df %>%
      arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
      group_by(run_id, social_distance_h) %>% 
      mutate(
        num_adopters_prop = num_adopters / N_nodes_actual,
        jump_in_adoption = num_adopters_prop - lag(num_adopters_prop, default = 0)
      ) %>%
      mutate(is_phase_transition_step = ifelse(jump_in_adoption >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)) %>%
      ungroup()
    
    first_transition_points <- phase_transition_data %>%
      filter(is_phase_transition_step == 1) %>%
      group_by(run_id, social_distance_h) %>%
      summarise(first_transition_IUL = min(innovation_iul_Gamma), .groups = 'drop')
    
    heatmap_data_list_temp <- list()
    for (iul_val in IUL_VALUES_SWEEP) {
      for (h_val in H_VALUES_SWEEP) {
        count_transitions <- first_transition_points %>%
          filter(social_distance_h == h_val, first_transition_IUL <= iul_val) %>% nrow()
        prop_transitions <- count_transitions / NUM_SEED_RUNS_TOTAL 
        heatmap_data_list_temp[[length(heatmap_data_list_temp) + 1]] <- data.frame(
          innovation_iul_Gamma = iul_val, social_distance_h = h_val,
          proportion_with_transition = prop_transitions,
          tau_mean_param = current_threshold_mean # Keep for potential internal faceting if needed
        )
      }
    }
    current_sd_all_means_heatmap_data[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- bind_rows(heatmap_data_list_temp)
    cat(paste0("    Phase transition analysis complete for (Mean τ = ", current_threshold_mean, ", SD τ = ", current_tau_sd, ").\n"))
    
  } # End inner loop for THRESHOLD_MEAN_SWEEP_LIST
  
  # --- Combine heatmap data for the CURRENT SD (from all its means) ---
  # And add the SD parameter for faceting the final plot
  if (length(current_sd_all_means_heatmap_data) > 0) {
    valid_heatmap_data_for_sd <- current_sd_all_means_heatmap_data[!sapply(current_sd_all_means_heatmap_data, is.null)]
    if (length(valid_heatmap_data_for_sd) > 0) {
      df_heatmap_this_sd <- bind_rows(valid_heatmap_data_for_sd)
      df_heatmap_this_sd$tau_sd_param <- current_tau_sd 
      all_sds_heatmap_data_list[[paste0("sd_", sprintf("%.2f", current_tau_sd))]] <- df_heatmap_this_sd
      
      # --- Generate and Save Plot for CURRENT SD ---
      df_heatmap_this_sd$social_distance_h_factor <- factor(sprintf("%.2f", df_heatmap_this_sd$social_distance_h))
      df_heatmap_this_sd$tau_mean_facet_label <- factor(paste0("Mean τ = ", sprintf("%.2f", df_heatmap_this_sd$tau_mean_param)))
      
      plot_for_current_sd <- ggplot(df_heatmap_this_sd, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = proportion_with_transition)) +
        geom_tile(color = "white", lwd = 0.15) +
        scale_fill_viridis_c(name = "Prop. of Runs\nwith Phase Transition\n(Jump >= 0.5)", limits = c(0, 1), option="plasma", n.breaks=6) +
        facet_wrap(~tau_mean_facet_label, ncol = 2) + 
        labs(
          x = expression(paste("Intrinsic Innovation Utility (", Gamma, ")")),
          y = "Maximum Social Distance (h)",
          title = paste("Phase Transition Maps for", CURRENT_GRAPH_TYPE_LABEL, "- Tau SD =", sprintf("%.2f", current_tau_sd)),
          subtitle = paste("Thresholds ~ N(μ=variable, σ=", sprintf("%.2f", current_tau_sd), "), ", NUM_SEED_RUNS_TOTAL, " runs per (Γ,h) per panel", sep="")
        ) +
        scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0,0)) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size=12),
          plot.subtitle = element_text(hjust = 0.5, size=8),
          axis.text.x = element_text(angle = 45, hjust = 1, size=7),
          axis.text.y = element_text(size=7), legend.position = "right",
          strip.text = element_text(face="bold", size=9) 
        )
      
      plot_filename_current_sd <- paste0(PLOTS_DIR, "phase_transition_sd", sprintf("%.2f", current_tau_sd), "_means03-06_", SEEDING_STRATEGY_FIXED, ".pdf")
      ggsave(plot_filename_current_sd, plot_for_current_sd, width = 10, height = 8) # Adjusted for potentially 4 means
      cat(paste0("    Saved plot: ", plot_filename_current_sd, "\n"))
    }
  }
  
  # --- Save Raw Data for CURRENT SD ---
  raw_data_filename_current_sd <- paste0(RESULTS_DIR, "phase_transition_raw_sd", sprintf("%.2f", current_tau_sd), "_means03-06_", SEEDING_STRATEGY_FIXED, ".rds")
  saveRDS(current_sd_all_means_raw_results, raw_data_filename_current_sd) # Save list of DFs (one per mean)
  cat(paste0("    Saved raw data: ", raw_data_filename_current_sd, "\n"))
  
  # Store the raw results for this SD in the grand list
  all_sds_raw_results_list[[paste0("sd_", sprintf("%.2f", current_tau_sd))]] <- current_sd_all_means_raw_results
  
} # End outer loop for TAU_NORMAL_SD_SWEEP_LIST
time_fin <- Sys.time()
time_total_parallel_3 <- difftime(time_fin, time_init, units = "auto")

cat("\nAll simulation sweeps complete.\n")

# -----------------------------------------------------------------------------
# 6. Final Grand Data Saving (optional, as data is saved per SD)
# -----------------------------------------------------------------------------
cat("Saving grand combined data objects...\n")

# Grand combined heatmap data (all SDs and all Means)
grand_combined_heatmap_df <- bind_rows(all_sds_heatmap_data_list)
saveRDS(grand_combined_heatmap_df, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_heatmap_data_all_sds_means_", SEEDING_STRATEGY_FIXED, ".rds"))

# Grand combined raw results (list of lists)
saveRDS(all_sds_raw_results_list, paste0(RESULTS_DIR, "phase_transition_GRAND_COMBINED_raw_results_all_sds_means_", SEEDING_STRATEGY_FIXED, ".rds"))

cat("All data saved. Script finished.\n")
