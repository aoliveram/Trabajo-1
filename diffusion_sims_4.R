# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Script to simulate complex contagion, identify phase transitions,
#              and visualize them as heatmaps of IUL vs. h for different
#              threshold distribution MEANS (fixed SD).
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

CURRENT_GRAPH_TYPE_LABEL <- "ATP-net"
NETWORKS_DIR <- "trabajo_1_files/ATP_network_ergm/"
NUM_NETWORK_INSTANCES_TO_USE <- 32 
N_NODES_GLOBAL <- 1000

# Fixed Standard Deviation for Threshold Distribution (τ_i)
TAU_NORMAL_SD_FIXED <- 0.16 # Punto 1: SD fija

# Means for Threshold Distribution (τ_i) - NOW A LIST TO SWEEP
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6) # Punto 1: Medias a variar

SEEDING_STRATEGY_FIXED <- "random"
NUM_SEED_RUNS_TOTAL <- NUM_NETWORK_INSTANCES_TO_USE # Cada run usa una red diferente

PHASE_TRANSITION_THRESHOLD_JUMP <- 0.50
NUM_CORES_TO_USE <- 8

# -----------------------------------------------------------------------------
# 2. Parameters to Sweep (IUL and h)
# -----------------------------------------------------------------------------
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP <- seq(0/12, 12/12, by = 1/12) # Punto 2: Granularidad de h

# -----------------------------------------------------------------------------
# 3. Main Simulation Loop
# -----------------------------------------------------------------------------
cat("Starting phase transition simulations...\n")

all_runs_results_list <- list() # Will store dataframes from each MEAN_TAU run

# Outer loop for different MEANS of Tau
total_mean_iterations <- length(THRESHOLD_MEAN_SWEEP_LIST)
for (mean_idx in 1:total_mean_iterations) {
  current_threshold_mean <- THRESHOLD_MEAN_SWEEP_LIST[mean_idx]
  
  cat(paste0("\n====================================================================\n"))
  cat(paste0("Processing for Tau Distribution: Normal(μ=", current_threshold_mean, ", σ=", TAU_NORMAL_SD_FIXED, ")\n"))
  cat(paste0(" (Iteration ", mean_idx, " of ", total_mean_iterations, " for threshold means)\n"))
  cat(paste0("====================================================================\n"))
  
  cl <- makeCluster(NUM_CORES_TO_USE, type = "FORK")
  registerDoParallel(cl)
  cat(paste0("  Registered ", NUM_CORES_TO_USE, " parallel workers for this mean τ iteration.\n"))
  cat(paste0("  Starting ", NUM_SEED_RUNS_TOTAL, " simulation runs (network+seed pairs)...\n"))
  
  list_of_results_for_this_mean <- foreach(
    run_idx = 1:NUM_SEED_RUNS_TOTAL,
    .combine = 'list',
    .multicombine = TRUE,
    .packages = c('igraph', 'dplyr', 'readr', 'intergraph', 'cluster'),
    .export = c('sweep_homoph_parameter', 'get_complex_plot',
                'NETWORKS_DIR', 'current_threshold_mean', 'TAU_NORMAL_SD_FIXED', # Pass current_threshold_mean
                'IUL_VALUES_SWEEP', 'H_VALUES_SWEEP', 'SEEDING_STRATEGY_FIXED',
                'N_NODES_GLOBAL' # Define this if used or remove if N_NODES_SPECIFIC_GRAPH is always used
    ), 
    .errorhandling = 'pass'
  ) %dopar% {
    
    worker_message_prefix <- paste0("    Worker ", Sys.getpid(), " [Mean τ:", sprintf("%.2f", current_threshold_mean), "] (Run ", run_idx, "/", NUM_SEED_RUNS_TOTAL, "): ")
    # cat(paste0(worker_message_prefix, "Starting.\n")) # Can be too verbose
    
    network_file_idx <- run_idx 
    current_network_path <- paste0(NETWORKS_DIR, "ATP_network_simulated_1000_mur_", sprintf("%03d", network_file_idx), ".rds")
    
    if (!file.exists(current_network_path)) {
      cat(paste0(worker_message_prefix, "ERROR - Network file not found: ", basename(current_network_path), ". Skipping.\n"))
      return(NULL) 
    }
    # cat(paste0(worker_message_prefix, "Loading network ", basename(current_network_path), "...\n")) # Verbose
    graph_for_this_run_ergm <- readRDS(current_network_path)
    graph_for_this_run <- asIgraph(graph_for_this_run_ergm)
    N_NODES_SPECIFIC_GRAPH <- vcount(graph_for_this_run)
    
    node_mur_q_specific <- V(graph_for_this_run)$q_i
    node_degrees_specific <- igraph::degree(graph_for_this_run)
    
    attributes_for_distance_specific <- data.frame(
      age = V(graph_for_this_run)$age, educ_num = V(graph_for_this_run)$educ_num,
      race = as.factor(V(graph_for_this_run)$race), relig = as.factor(V(graph_for_this_run)$relig),
      sex = as.factor(V(graph_for_this_run)$sex)
    )
    social_distance_matrix_d_ij_specific <- as.matrix(daisy(attributes_for_distance_specific, metric = "gower"))
    
    set.seed(run_idx * 1000 + round(current_threshold_mean * 100)) 
    node_thresholds_tau_frac_specific <- rnorm(
      n = N_NODES_SPECIFIC_GRAPH, mean = current_threshold_mean, sd = TAU_NORMAL_SD_FIXED
    )
    node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific < 0] <- 0
    node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific > 1] <- 1
    
    node_thresholds_count_for_cluster_specific <- round(node_thresholds_tau_frac_specific * node_degrees_specific)
    node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific > 0] <- 1
    node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific == 0] <- 0
    
    set.seed(run_idx * 2000) 
    primary_seed_for_this_run <- as.integer(sample(V(graph_for_this_run), 1))
    
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
    # cat(paste0(worker_message_prefix, "Seed node: ", primary_seed_for_this_run, ", Initial cluster size: ", length(initial_infectors_for_this_sim_run), "\n")) # Verbose
    
    # cat(paste0(worker_message_prefix, "Sweeping IUL and h parameters...\n")) # Verbose
    df_one_full_run <- sweep_homoph_parameter(
      primary_seed_id_arg = primary_seed_for_this_run, 
      N_nodes_arg = N_NODES_SPECIFIC_GRAPH,
      graph_obj_arg = graph_for_this_run,
      node_individual_thresholds_tau_arg = node_thresholds_tau_frac_specific,
      node_mur_q_arg = node_mur_q_specific,
      all_innovation_iul_Gamma_values = IUL_VALUES_SWEEP,
      all_social_distance_h_values = H_VALUES_SWEEP,
      initial_infectors_vector_arg = initial_infectors_for_this_sim_run,
      d_ij_matrix = social_distance_matrix_d_ij_specific
    )
    
    df_one_full_run$run_id <- run_idx 
    df_one_full_run$network_instance_file <- basename(current_network_path)
    df_one_full_run$threshold_mean_param <- current_threshold_mean # The mean used for this run
    df_one_full_run$threshold_sd_param <- TAU_NORMAL_SD_FIXED 
    df_one_full_run$N_nodes_actual <- N_NODES_SPECIFIC_GRAPH
    
    # Simple progress indication from worker (optional, can be too much output)
    # if (run_idx %% (NUM_SEED_RUNS_TOTAL / 10) == 0 || run_idx == NUM_SEED_RUNS_TOTAL) { # ~ every 10%
    #   cat(paste0(worker_message_prefix, "Progress mark. Rows generated: ", nrow(df_one_full_run), "\n"))
    # }
    return(df_one_full_run)
  }
  
  stopCluster(cl)
  cat(paste0("  Parallel simulations for Mean τ = ", current_threshold_mean, " finished (", length(list_of_results_for_this_mean) , " runs collected).\n"))
  
  current_mean_results_raw_df_unfiltered <- bind_rows(list_of_results_for_this_mean)
  # Filter out NULLs more robustly
  valid_results_indices <- !sapply(list_of_results_for_this_mean, function(x) inherits(x, "simpleError") || is.null(x))
  if(sum(valid_results_indices) > 0) {
    current_mean_results_raw_df <- bind_rows(list_of_results_for_this_mean[valid_results_indices])
  } else {
    current_mean_results_raw_df <- data.frame() # Empty dataframe
  }
  
  if (nrow(current_mean_results_raw_df) == 0) {
    cat(paste0("  WARNING: No valid results obtained for Mean τ = ", current_threshold_mean, ". Skipping analysis for this mean.\n"))
    next 
  }
  
  all_runs_results_list[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- current_mean_results_raw_df
  cat(paste0("  Stored results for Mean τ = ", current_threshold_mean, ". Total rows: ", nrow(current_mean_results_raw_df),"\n"))
  
} 

cat("\nAll simulation iterations complete.\n")

# -----------------------------------------------------------------------------
# 4. Phase Transition Identification & Data Aggregation for Plotting
# -----------------------------------------------------------------------------
cat("Identifying phase transitions across all collected results...\n")

list_of_heatmaps_dfs <- list()

for (mean_label in names(all_runs_results_list)) {
  current_mean_raw_df <- all_runs_results_list[[mean_label]]
  current_mean_value <- as.numeric(gsub("mean_", "", mean_label)) 
  
  cat(paste0("  Processing phase transitions for Mean τ: ", current_mean_value, "\n"))
  
  if (is.null(current_mean_raw_df) || nrow(current_mean_raw_df) == 0) {
    cat(paste0("    No data for Mean τ = ", current_mean_value, ". Skipping heatmap generation.\n"))
    next
  }
  
  phase_transition_data_this_mean <- current_mean_raw_df %>%
    arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
    group_by(run_id, social_distance_h) %>% 
    mutate(
      num_adopters_prop = num_adopters / N_nodes_actual,
      jump_in_adoption = num_adopters_prop - lag(num_adopters_prop, default = 0)
    ) %>%
    mutate(
      is_phase_transition_step = ifelse(jump_in_adoption >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
    ) %>%
    ungroup()
  
  first_transition_points_this_mean <- phase_transition_data_this_mean %>%
    filter(is_phase_transition_step == 1) %>%
    group_by(run_id, social_distance_h) %>%
    summarise(
      first_transition_IUL = min(innovation_iul_Gamma),
      .groups = 'drop'
    )
  
  heatmap_data_list_this_mean <- list()
  for (iul_val in IUL_VALUES_SWEEP) {
    for (h_val in H_VALUES_SWEEP) {
      count_transitions_at_or_before_IUL <- first_transition_points_this_mean %>%
        filter(social_distance_h == h_val, first_transition_IUL <= iul_val) %>%
        nrow() 
      
      prop_transitions <- count_transitions_at_or_before_IUL / NUM_SEED_RUNS_TOTAL 
      
      heatmap_data_list_this_mean[[length(heatmap_data_list_this_mean) + 1]] <- data.frame(
        innovation_iul_Gamma = iul_val,
        social_distance_h = h_val,
        proportion_with_transition = prop_transitions,
        tau_mean_param = current_mean_value 
      )
    }
  }
  list_of_heatmaps_dfs[[mean_label]] <- bind_rows(heatmap_data_list_this_mean)
  cat(paste0("    Finished processing for Mean τ: ", current_mean_value, "\n"))
}

combined_heatmap_df <- bind_rows(list_of_heatmaps_dfs)

# -----------------------------------------------------------------------------
# 5. Visualization (Faceted Heatmap by Mean Tau)
# -----------------------------------------------------------------------------
cat("Generating faceted phase transition heatmap...\n")

if (nrow(combined_heatmap_df) > 0) {
  combined_heatmap_df$social_distance_h_factor <- factor(sprintf("%.2f", combined_heatmap_df$social_distance_h))
  combined_heatmap_df$tau_mean_facet_label <- factor(paste0("Mean τ = ", sprintf("%.2f", combined_heatmap_df$tau_mean_param)))
  
  faceted_phase_transition_plot <- ggplot(combined_heatmap_df, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = proportion_with_transition)) +
    geom_tile(color = "white", lwd = 0.15) +
    scale_fill_viridis_c(name = "Prop. of Runs\nwith Phase Transition\n(Jump >= 0.5)", limits = c(0, 1), option="plasma", n.breaks=6) +
    facet_wrap(~tau_mean_facet_label, ncol = 2) + 
    labs(
      x = expression(paste("Intrinsic Innovation Utility (", Gamma, ")")),
      y = "Maximum Social Distance (h)",
      title = paste("Phase Transition Maps for", CURRENT_GRAPH_TYPE_LABEL),
      subtitle = paste("Thresholds ~ N(μ=variable, σ=", TAU_NORMAL_SD_FIXED, "), ", NUM_SEED_RUNS_TOTAL, " runs (network+seed) per (Γ,h) point per panel", sep="")
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0,0)) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size=12),
      plot.subtitle = element_text(hjust = 0.5, size=8),
      axis.text.x = element_text(angle = 45, hjust = 1, size=7),
      axis.text.y = element_text(size=7),
      legend.position = "right",
      strip.text = element_text(face="bold", size=9) 
    )
  
  print(faceted_phase_transition_plot)
} else {
  cat("No data available to generate faceted heatmap.\n")
}

# Save the combined raw data and combined heatmap data
saveRDS(all_runs_results_list, paste0("trabajo_1_files/diffusion_simulation_files/phase_transition_means030-060_sd",TAU_NORMAL_SD_FIXED,".rds"))
ggsave(paste0("trabajo_1_plots/diffusion_simulation_plots/phase_transition_means030-060_sd",TAU_NORMAL_SD_FIXED, ".pdf"),
       faceted_phase_transition_plot, 
       width = 8, height = 5.5)
