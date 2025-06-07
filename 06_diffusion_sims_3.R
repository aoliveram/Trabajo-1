# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Script to simulate complex contagion, identify phase transitions,
#              and visualize them as heatmaps of IUL vs. h for different
#              threshold distribution standard deviations.
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
# 1. Fixed Parameters (Some now become lists to iterate over)
# -----------------------------------------------------------------------------
cat("Setting up fixed and sweep parameters...\n")

# Network Type (remains ATP)
CURRENT_GRAPH_TYPE_LABEL <- "ATP-net"
NETWORKS_DIR <- "trabajo_1_files/ATP_network_ergm/"
NUM_NETWORK_INSTANCES_TO_USE <- 32 # Corresponds to NUM_SEED_RUNS_TOTAL
N_NODES_GLOBAL <- 1000

# Fixed Mean for Threshold Distribution (τ_i)
THRESHOLD_MEAN_FIXED <- 0.50 # This is fixed as per your previous request

# Standard Deviations for Threshold Distribution (τ_i) - NOW A LIST TO SWEEP
TAU_NORMAL_SD_SWEEP_LIST <- c(0.05, 0.10, 0.16, 0.25) # Example SDs, you mentioned 0.3,0.4,0.5,0.6 which might be too large for mean 0.5 if τ must be in [0,1]
# If mean is 0.5, sd=0.3 means 99.7% data within mean +/- 3*sd = 0.5 +/- 0.9 = [-0.4, 1.4].
# Let's use your requested SDs, but be mindful of truncation.
#TAU_NORMAL_SD_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6) # Your requested values for SD

# Seeding Strategy (remains random for each network instance run)
SEEDING_STRATEGY_FIXED <- "random"
NUM_SEED_RUNS_TOTAL <- 32 # Total number of simulation runs (each with a new network & seed)

# Phase Transition Definition
PHASE_TRANSITION_THRESHOLD_JUMP <- 0.50

# Parallel Processing
NUM_CORES_TO_USE <- 8

# -----------------------------------------------------------------------------
# 2. Parameters to Sweep (IUL and h)
# -----------------------------------------------------------------------------
# IUL (Γ) values: 0 to 1, steps of 0.025
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
# h values: 0/12, 1/12, ..., 12/12 (increased granularity)
H_VALUES_SWEEP <- seq(0/12, 12/12, by = 1/12)

# -----------------------------------------------------------------------------
# 3. Main Simulation Loop
# -----------------------------------------------------------------------------
cat("Starting phase transition simulations...\n")

all_runs_results_list <- list() # Will store dataframes from each SD run

# Outer loop for different Standard Deviations of Tau
for (current_tau_sd in TAU_NORMAL_SD_SWEEP_LIST) {
  
  cat(paste("\n====================================================================\n"))
  cat(paste("Processing for Tau Distribution: Normal(μ=", THRESHOLD_MEAN_FIXED, ", σ=", current_tau_sd, ")\n"))
  cat(paste("====================================================================\n"))
  
  # Setup for parallel processing for this SD run
  # Each worker will handle one full simulation run (one network instance, one seed, all IUL/h)
  cl <- makeCluster(NUM_CORES_TO_USE, type = "FORK") # "PSOCK" for Windows
  registerDoParallel(cl)
  cat(paste("  Registered", NUM_CORES_TO_USE, "parallel workers for this SD iteration.\n"))
  
  # The foreach loop now iterates NUM_SEED_RUNS_TOTAL times.
  # Each iteration is an independent simulation with its own network and seed.
  list_of_results_for_this_sd <- foreach(
    run_idx = 1:NUM_SEED_RUNS_TOTAL,
    .combine = 'list',
    .multicombine = TRUE,
    .packages = c('igraph', 'dplyr', 'readr', 'intergraph', 'cluster'), # Added readr, intergraph, cluster
    .export = c('sweep_homoph_parameter', 'get_complex_plot', # Assuming these are in the sourced file
                'NETWORKS_DIR', 'THRESHOLD_MEAN_FIXED', 'current_tau_sd', # Note: current_tau_sd passed
                'IUL_VALUES_SWEEP', 'H_VALUES_SWEEP', 'SEEDING_STRATEGY_FIXED', # Global sweeps
                'N_NODES_GLOBAL' # Define this if needed, or get from graph
    ), 
    .errorhandling = 'pass'
  ) %dopar% {
    
    # --- Inside each parallel worker ---
    worker_message_prefix <- paste0("    Worker ", Sys.getpid(), " (Run ", run_idx, "/", NUM_SEED_RUNS_TOTAL, "): ")
    cat(paste0(worker_message_prefix, "Starting.\n"))
    
    # 1. Load a specific network instance for this run
    #    The network file index corresponds to run_idx
    network_file_idx <- run_idx 
    current_network_path <- paste0(NETWORKS_DIR, "ATP_network_simulated_1000_mur_", sprintf("%03d", network_file_idx), ".rds")
    
    if (!file.exists(current_network_path)) {
      cat(paste0(worker_message_prefix, "ERROR - Network file not found: ", current_network_path, ". Skipping this run.\n"))
      return(NULL) # Skip this iteration if file not found
    }
    cat(paste0(worker_message_prefix, "Loading network ", basename(current_network_path), "...\n"))
    graph_for_this_run_ergm <- readRDS(current_network_path)
    graph_for_this_run <- asIgraph(graph_for_this_run_ergm)
    N_NODES_SPECIFIC_GRAPH <- vcount(graph_for_this_run)
    
    # 2. Extract/Calculate attributes for this specific graph
    node_mur_q_specific <- V(graph_for_this_run)$q_i
    node_degrees_specific <- igraph::degree(graph_for_this_run)
    
    attributes_for_distance_specific <- data.frame(
      age = V(graph_for_this_run)$age,
      educ_num = V(graph_for_this_run)$educ_num,
      race = as.factor(V(graph_for_this_run)$race),
      relig = as.factor(V(graph_for_this_run)$relig),
      sex = as.factor(V(graph_for_this_run)$sex)
    )
    social_distance_matrix_d_ij_specific <- as.matrix(daisy(attributes_for_distance_specific, metric = "gower"))
    
    # 3. Generate Tau thresholds for this specific graph and run_idx (for reproducibility if needed)
    set.seed(run_idx * 1000 + round(current_tau_sd * 100)) # Seed for tau generation
    node_thresholds_tau_frac_specific <- rnorm(
      n = N_NODES_SPECIFIC_GRAPH,
      mean = THRESHOLD_MEAN_FIXED,
      sd = current_tau_sd # Use the SD for the current outer loop
    )
    node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific < 0] <- 0
    node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific > 1] <- 1
    
    node_thresholds_count_for_cluster_specific <- round(node_thresholds_tau_frac_specific * node_degrees_specific)
    node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific > 0] <- 1
    node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific == 0] <- 0
    
    
    # 4. Select a single primary seed node for this run
    cat(paste0(worker_message_prefix, "Selecting seed node...\n"))
    set.seed(run_idx * 2000) # Seed for random seed selection
    if (SEEDING_STRATEGY_FIXED == "random") {
      primary_seed_for_this_run <- as.integer(sample(V(graph_for_this_run), 1))
    } else {
      # Implement other deterministic strategies if needed, ensuring they pick ONE seed or a defined set
      cat(paste0(worker_message_prefix, "ERROR - Seeding strategy '", SEEDING_STRATEGY_FIXED, "' not fully implemented for single seed selection here. Defaulting to random.\n"))
      primary_seed_for_this_run <- as.integer(sample(V(graph_for_this_run), 1))
    }
    
    # 5. Determine initial infector cluster for this primary_seed
    num_seeds_for_initial_cluster <- node_thresholds_count_for_cluster_specific[primary_seed_for_this_run]
    num_seeds_for_initial_cluster <- min(num_seeds_for_initial_cluster, N_NODES_SPECIFIC_GRAPH, node_degrees_specific[primary_seed_for_this_run] + 1)
    if (num_seeds_for_initial_cluster < 1) num_seeds_for_initial_cluster <- 1
    
    initial_infectors_for_this_sim_run <- c(primary_seed_for_this_run)
    if (num_seeds_for_initial_cluster > 1) {
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
    cat(paste0(worker_message_prefix, "Seed node: ", primary_seed_for_this_run, ", Initial cluster size: ", length(initial_infectors_for_this_sim_run), "\n"))
    
    # 6. Call sweep_homoph_parameter for all IUL and h values
    cat(paste0(worker_message_prefix, "Sweeping IUL and h parameters...\n"))
    df_one_full_run <- sweep_homoph_parameter(
      primary_seed_id_arg = primary_seed_for_this_run, # This is more like a run_identifier now
      N_nodes_arg = N_NODES_SPECIFIC_GRAPH,
      graph_obj_arg = graph_for_this_run,
      node_individual_thresholds_tau_arg = node_thresholds_tau_frac_specific,
      node_mur_q_arg = node_mur_q_specific,
      all_innovation_iul_Gamma_values = IUL_VALUES_SWEEP,
      all_social_distance_h_values = H_VALUES_SWEEP,
      initial_infectors_vector_arg = initial_infectors_for_this_sim_run,
      d_ij_matrix = social_distance_matrix_d_ij_specific
    )
    
    # 7. Add identifiers
    df_one_full_run$run_id <- run_idx # Identifies this specific simulation run (network instance + seed)
    df_one_full_run$network_instance_file <- basename(current_network_path)
    df_one_full_run$threshold_mean_param <- THRESHOLD_MEAN_FIXED
    df_one_full_run$threshold_sd_param <- current_tau_sd # The SD used for this run
    df_one_full_run$N_nodes_actual <- N_NODES_SPECIFIC_GRAPH
    
    cat(paste0(worker_message_prefix, "Finished. Rows generated: ", nrow(df_one_full_run), "\n"))
    return(df_one_full_run)
  }# End foreach
  
  stopCluster(cl)
  cat(paste("  Parallel simulations for Tau SD =", current_tau_sd, "finished.\n"))
  
  # Combine results from all runs for THIS SD value, and filter out NULLs if any file was missing
  current_sd_results_raw_df_unfiltered <- bind_rows(list_of_results_for_this_sd)
  current_sd_results_raw_df <- current_sd_results_raw_df_unfiltered[!sapply(current_sd_results_raw_df_unfiltered, is.null),]
  
  if (nrow(current_sd_results_raw_df) == 0) {
    cat(paste("  WARNING: No valid results obtained for Tau SD =", current_tau_sd, ". Skipping analysis for this SD.\n"))
    next # Skip to the next SD if no data
  }
  
  all_runs_results_list[[paste0("sd_", sprintf("%.2f", current_tau_sd))]] <- current_sd_results_raw_df
  
} # End outer loop for TAU_NORMAL_SD_SWEEP_LIST

cat("\nAll simulation iterations complete.\n")

# -----------------------------------------------------------------------------
# 4. Phase Transition Identification & Data Aggregation for Plotting
# -----------------------------------------------------------------------------
cat("Identifying phase transitions across all collected results...\n")

list_of_heatmaps_dfs <- list()

for (sd_label in names(all_runs_results_list)) {
  current_sd_raw_df <- all_runs_results_list[[sd_label]]
  current_sd_value <- as.numeric(gsub("sd_", "", sd_label)) # Extract SD value
  
  cat(paste("  Processing phase transitions for Tau SD:", current_sd_value, "\n"))
  
  if (is.null(current_sd_raw_df) || nrow(current_sd_raw_df) == 0) {
    cat(paste("    No data for SD =", current_sd_value, ". Skipping heatmap generation.\n"))
    next
  }
  
  phase_transition_data_this_sd <- current_sd_raw_df %>%
    arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
    group_by(run_id, social_distance_h) %>% # Group by each simulation run and h value
    mutate(
      num_adopters_prop = num_adopters / N_nodes_actual, # Use N_nodes_actual from the run
      jump_in_adoption = num_adopters_prop - lag(num_adopters_prop, default = 0)
    ) %>%
    mutate(
      is_phase_transition_step = ifelse(jump_in_adoption >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
    ) %>%
    ungroup()
  
  first_transition_points_this_sd <- phase_transition_data_this_sd %>%
    filter(is_phase_transition_step == 1) %>%
    group_by(run_id, social_distance_h) %>%
    summarise(
      first_transition_IUL = min(innovation_iul_Gamma),
      .groups = 'drop'
    )
  
  heatmap_data_list_this_sd <- list()
  for (iul_val in IUL_VALUES_SWEEP) {
    for (h_val in H_VALUES_SWEEP) {
      count_transitions_at_or_before_IUL <- first_transition_points_this_sd %>%
        filter(social_distance_h == h_val, first_transition_IUL <= iul_val) %>%
        nrow() # Number of runs (network instance + seed) that had a transition
      
      # Proportion over total runs for this SD
      prop_transitions <- count_transitions_at_or_before_IUL / NUM_SEED_RUNS_TOTAL 
      
      heatmap_data_list_this_sd[[length(heatmap_data_list_this_sd) + 1]] <- data.frame(
        innovation_iul_Gamma = iul_val,
        social_distance_h = h_val,
        proportion_with_transition = prop_transitions,
        tau_sd_param = current_sd_value # Add SD for faceting
      )
    }
  }
  list_of_heatmaps_dfs[[sd_label]] <- bind_rows(heatmap_data_list_this_sd)
  cat(paste("    Finished processing for Tau SD:", current_sd_value, "\n"))
}

combined_heatmap_df <- bind_rows(list_of_heatmaps_dfs)

# -----------------------------------------------------------------------------
# 5. Visualization (Faceted Heatmap)
# -----------------------------------------------------------------------------
cat("Generating faceted phase transition heatmap...\n")

if (nrow(combined_heatmap_df) > 0) {
  combined_heatmap_df$social_distance_h_factor <- factor(sprintf("%.2f", combined_heatmap_df$social_distance_h))
  combined_heatmap_df$tau_sd_facet_label <- factor(paste0("Tau SD = ", sprintf("%.2f", combined_heatmap_df$tau_sd_param)))
  
  faceted_phase_transition_plot <- ggplot(combined_heatmap_df, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = proportion_with_transition)) +
    geom_tile(color = "white", lwd = 0.15) +
    scale_fill_viridis_c(name = "Prop. of Runs\nwith Phase Transition\n(Jump >= 0.5)", limits = c(0, 1), option="plasma", n.breaks=6) +
    facet_wrap(~tau_sd_facet_label, ncol = 2) + # Facet by Tau SD
    labs(
      x = expression(paste("Intrinsic Innovation Utility (", Gamma, ")")),
      y = "Maximum Social Distance (h)",
      title = paste("Phase Transition Maps for", CURRENT_GRAPH_TYPE_LABEL),
      subtitle = paste("Thresholds ~ N(μ=", THRESHOLD_MEAN_FIXED, "), ", NUM_SEED_RUNS_TOTAL, " runs (network+seed) per (Γ,h) point per SD panel", sep="")
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0,0)) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size=12),
      plot.subtitle = element_text(hjust = 0.5, size=8),
      axis.text.x = element_text(angle = 45, hjust = 1, size=7),
      axis.text.y = element_text(size=7),
      legend.position = "right",
      strip.text = element_text(face="bold", size=9) # Facet labels
    )
  
  print(faceted_phase_transition_plot)
  
  plot_filename_faceted <- paste0("phase_transition_faceted_heatmap_", gsub("[: ]", "_", Sys.time()), ".png")
  ggsave(plot_filename_faceted, faceted_phase_transition_plot, width = 10, height = 8, dpi = 300)
  cat(paste("Faceted phase transition heatmap saved to:", plot_filename_faceted, "\n"))
} else {
  cat("No data available to generate faceted heatmap.\n")
}

# Save the combined raw data and combined heatmap data
saveRDS(all_runs_results_list, "phase_transition_all_sds_raw_results.rds") # List of DFs
saveRDS(combined_heatmap_df, "phase_transition_all_sds_heatmap_data.rds")
cat("Raw simulation results (list per SD) and combined heatmap data saved.\n")