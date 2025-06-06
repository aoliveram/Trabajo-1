# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Script to simulate complex contagion, identify phase transitions,
#              and visualize them as a heatmap of IUL vs. h.
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
# 1. Fixed Parameters and Network Loading
# -----------------------------------------------------------------------------
cat("Setting up fixed parameters and loading network...\n")

# Network Type
CURRENT_GRAPH_TYPE_LABEL <- "ATP-net"
ATP_NETWORK <- TRUE # For loading ATP networks

# Load ATP Network (using only the first instance for this focused study)
# If you want to average over multiple ATP-ERGM instances, you'd loop this or load all
# and then loop the main simulation over them. For now, one instance.
if (ATP_NETWORK) {
  networks_dir <- "trabajo_1_files/ATP_network_ergm/"
  # Assuming you want to run on one representative ATP network for this specific phase transition study
  # If you have 5 and want to run on all, you'll need an outer loop for graphs_list_to_simulate
  atp_net_path <- paste0(networks_dir, "ATP_network_simulated_1000_mur_001.rds") # Using instance 1
  if (!file.exists(atp_net_path)) {
    stop(paste("ATP network file not found:", atp_net_path))
  }
  single_atp_graph <- readRDS(atp_net_path)
  single_atp_graph <- asIgraph(single_atp_graph)
  graphs_list_to_simulate <- list(single_atp_graph) # List with one graph
  N_NODES_GLOBAL <- vcount(graphs_list_to_simulate[[1]])
} else {
  stop("This script is configured for ATP_NETWORK = TRUE only for phase transition analysis.")
}

# Fixed Node Attributes and Social Distance
# These are calculated once per graph instance
current_graph_obj_fixed <- graphs_list_to_simulate[[1]] # Only one graph in the list now
NODE_MUR_Q_FIXED <- V(current_graph_obj_fixed)$q_i
NODE_DEGREES_FIXED <- igraph::degree(current_graph_obj_fixed)

attributes_for_distance_fixed <- data.frame(
  age = V(current_graph_obj_fixed)$age,
  educ_num = V(current_graph_obj_fixed)$educ_num,
  race = as.factor(V(current_graph_obj_fixed)$race),
  relig = as.factor(V(current_graph_obj_fixed)$relig),
  sex = as.factor(V(current_graph_obj_fixed)$sex)
)
SOCIAL_DISTANCE_MATRIX_D_IJ_FIXED <- as.matrix(daisy(attributes_for_distance_fixed, metric = "gower"))

# Fixed Threshold Distribution (τ_i)
THRESHOLD_MEAN_FIXED <- 0.50 # As per your point 5.3
TAU_NORMAL_DISTRIBUTION_SD_FIXED <- 0.16 # As per your point 5.3
set.seed(42) # For reproducibility of threshold assignment
NODE_THRESHOLDS_TAU_FRAC_FIXED <- rnorm(
  n = N_NODES_GLOBAL,
  mean = THRESHOLD_MEAN_FIXED,
  sd = TAU_NORMAL_DISTRIBUTION_SD_FIXED
)
NODE_THRESHOLDS_TAU_FRAC_FIXED[NODE_THRESHOLDS_TAU_FRAC_FIXED < 0] <- 0
NODE_THRESHOLDS_TAU_FRAC_FIXED[NODE_THRESHOLDS_TAU_FRAC_FIXED > 1] <- 1

# Fixed Count Thresholds (for initial cluster in seeding)
NODE_THRESHOLDS_COUNT_FOR_CLUSTER_FIXED <- round(NODE_THRESHOLDS_TAU_FRAC_FIXED * NODE_DEGREES_FIXED)
NODE_THRESHOLDS_COUNT_FOR_CLUSTER_FIXED[NODE_THRESHOLDS_COUNT_FOR_CLUSTER_FIXED <= 0 & NODE_THRESHOLDS_TAU_FRAC_FIXED > 0] <- 1
NODE_THRESHOLDS_COUNT_FOR_CLUSTER_FIXED[NODE_THRESHOLDS_COUNT_FOR_CLUSTER_FIXED <= 0 & NODE_THRESHOLDS_TAU_FRAC_FIXED == 0] <- 0


# Fixed Seeding Strategy Parameters
SEEDING_STRATEGY_FIXED <- "random" 
NUM_SEED_NODES_TO_TEST_FIXED <- 16 # Your point 5.1

# Phase Transition Definition
PHASE_TRANSITION_THRESHOLD_JUMP <- 0.50 # Your point 3 (e.g., jump >= 50% adopters)

# -----------------------------------------------------------------------------
# 2. Parameters to Sweep for Phase Transition Analysis
# -----------------------------------------------------------------------------
cat("Defining sweep parameters...\n")
# IUL (Γ) values: 0 to 1, steps of 0.025 (point 4.1)
IUL_values_sweep <- seq(0.0, 1.0, by = 0.025)
# h values: (0/6, ..., 6/6) (point 4.2)
h_values_sweep <- seq(0/6, 6/6, by = 1/6)

# -----------------------------------------------------------------------------
# 3. Main Simulation Loop for Phase Transition Analysis
# -----------------------------------------------------------------------------
cat("Starting phase transition simulations...\n")

# Select primary seed nodes ONCE based on the fixed strategy
if (SEEDING_STRATEGY_FIXED == "random") {
  set.seed(101) # For reproducibility of random seed selection
  primary_seed_nodes_for_study <- as.integer(sample(V(current_graph_obj_fixed), NUM_SEED_NODES_TO_TEST_FIXED, replace = FALSE))
} else if (SEEDING_STRATEGY_FIXED == "PLci_top") { # Example if you wanted to use PLci
  adj_mat_for_plci_calc <- as.matrix(as_adjacency_matrix(current_graph_obj_fixed))
  plci_values <- sapply(1:N_NODES_GLOBAL, function(node_idx_for_plci) {
    get_PLci( # Assuming get_PLci uses NODE_THRESHOLDS_COUNT_FOR_CLUSTER_FIXED
      seed_node = node_idx_for_plci, N_nodes = N_NODES_GLOBAL, graph_obj = current_graph_obj_fixed,
      adj_matrix = adj_mat_for_plci_calc,
      node_thresholds = NODE_THRESHOLDS_COUNT_FOR_CLUSTER_FIXED, 
      num_initial_cluster_seeds_for_node = NODE_THRESHOLDS_COUNT_FOR_CLUSTER_FIXED[node_idx_for_plci],
      model_output_list_placeholder = NULL
    )[5]
  })
  primary_seed_nodes_for_study <- order(as.numeric(plci_values), decreasing = TRUE)[1:min(N_NODES_GLOBAL, NUM_SEED_NODES_TO_TEST_FIXED)]
} else {
  stop("Unsupported SEEDING_STRATEGY_FIXED.")
}
cat(paste("  Selected", length(primary_seed_nodes_for_study), "primary seed nodes for the study using strategy:", SEEDING_STRATEGY_FIXED, "\n"))


# Setup for parallel processing

cl <- makeCluster(8, type = "FORK") # "PSOCK" for Windows
registerDoParallel(cl)
cat(paste("  Registered", num_cores_to_use, "parallel workers.\n"))

# Each element in this list will be a dataframe from one primary_seed_node,
# containing results for all IUL and h combinations.
list_of_results_per_seed_node <- foreach(
  current_primary_seed_id = primary_seed_nodes_for_study,
  .combine = 'list',
  .multicombine = TRUE, # Important for list combination
  .packages = c('igraph', 'dplyr'),
  .export = c('sweep_homoph_parameter', 'get_complex_plot', # Assuming sweep_homoph_parameter calls get_complex_plot
              'N_NODES_GLOBAL', 'current_graph_obj_fixed',
              'NODE_THRESHOLDS_TAU_FRAC_FIXED', 'NODE_THRESHOLDS_COUNT_FOR_CLUSTER_FIXED',
              'NODE_MUR_Q_FIXED', 'NODE_DEGREES_FIXED',
              'IUL_values_sweep', 'h_values_sweep', 'SOCIAL_DISTANCE_MATRIX_D_IJ_FIXED'),
  .errorhandling = 'pass'
) %dopar% {
  
  # Determine the actual initial infector cluster for this primary_seed_id
  num_seeds_for_initial_cluster <- NODE_THRESHOLDS_COUNT_FOR_CLUSTER_FIXED[current_primary_seed_id]
  num_seeds_for_initial_cluster <- min(num_seeds_for_initial_cluster, N_NODES_GLOBAL, NODE_DEGREES_FIXED[current_primary_seed_id] + 1)
  if (num_seeds_for_initial_cluster < 1) num_seeds_for_initial_cluster <- 1
  
  initial_infectors_for_this_run <- c(current_primary_seed_id)
  if (num_seeds_for_initial_cluster > 1) {
    neighbors_of_primary <- as.numeric(neighbors(current_graph_obj_fixed, current_primary_seed_id, mode="total"))
    num_additional_needed <- num_seeds_for_initial_cluster - 1
    
    if (length(neighbors_of_primary) >= num_additional_needed) {
      initial_infectors_for_this_run <- c(initial_infectors_for_this_run, sample(neighbors_of_primary, num_additional_needed))
    } else {
      initial_infectors_for_this_run <- c(initial_infectors_for_this_run, neighbors_of_primary)
      still_needed_more <- num_seeds_for_initial_cluster - length(initial_infectors_for_this_run)
      if (still_needed_more > 0) {
        potential_others <- setdiff(1:N_NODES_GLOBAL, initial_infectors_for_this_run)
        if (length(potential_others) > 0) {
          initial_infectors_for_this_run <- c(initial_infectors_for_this_run, sample(potential_others, min(length(potential_others), still_needed_more)))
        }
      }
    }
  }
  initial_infectors_for_this_run <- unique(initial_infectors_for_this_run)
  
  # Call sweep_homoph_parameter which iterates over IUL_values_sweep and h_values_sweep
  # and calls get_complex_plot for each combination.
  df_one_seed_all_IUL_h <- sweep_homoph_parameter(
    primary_seed_id_arg = current_primary_seed_id,
    N_nodes_arg = N_NODES_GLOBAL,
    graph_obj_arg = current_graph_obj_fixed,
    node_individual_thresholds_tau_arg = NODE_THRESHOLDS_TAU_FRAC_FIXED,
    node_mur_q_arg = NODE_MUR_Q_FIXED,
    all_innovation_iul_Gamma_values = IUL_values_sweep,
    all_social_distance_h_values = h_values_sweep,
    initial_infectors_vector_arg = initial_infectors_for_this_run,
    d_ij_matrix = SOCIAL_DISTANCE_MATRIX_D_IJ_FIXED
  )
  
  # Add identifiers to the dataframe
  df_one_seed_all_IUL_h$primary_seed_node <- current_primary_seed_id
  df_one_seed_all_IUL_h$graph_type <- CURRENT_GRAPH_TYPE_LABEL # Only one graph type
  df_one_seed_all_IUL_h$threshold_mean <- THRESHOLD_MEAN_FIXED
  df_one_seed_all_IUL_h$threshold_sd <- TAU_NORMAL_DISTRIBUTION_SD_FIXED
  
  return(df_one_seed_all_IUL_h)
}

stopCluster(cl)
cat("  Parallel simulations finished.\n")

# Combine results from all primary seed nodes
all_sim_results_raw_df <- bind_rows(list_of_results_per_seed_node)

# -----------------------------------------------------------------------------
# 4. Phase Transition Identification
# -----------------------------------------------------------------------------
cat("Identifying phase transitions...\n")

# Ensure IUL_values_sweep is sorted for diff() to work correctly
# The output from sweep_homoph_parameter should already have innovation_iul_Gamma sorted
# for each seed and h value.

phase_transition_data <- all_sim_results_raw_df %>%
  arrange(primary_seed_node, social_distance_h, innovation_iul_Gamma) %>% # Ensure correct order
  group_by(primary_seed_node, social_distance_h) %>%
  mutate(
    num_adopters_prop = num_adopters / N_NODES_GLOBAL,
    jump_in_adoption = num_adopters_prop - lag(num_adopters_prop, default = 0) # Could also use diff() but lag is safer for first element
  ) %>%
  # Identify if a phase transition occurred for this (IUL, h) step FOR THIS SEED
  mutate(
    is_phase_transition_step = ifelse(jump_in_adoption >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
  ) %>%
  ungroup()

# Now, for each (IUL, h) pair, calculate the proportion of SEEDS that showed a transition AT or BEFORE this IUL
# This is tricky. A simpler approach is: for each (IUL,h) what proportion of seeds *experienced* a jump >= threshold *at that IUL step*?
# Or, for each (h), at what IUL did a transition occur for each seed?

# Let's find for each (seed, h) combination, the IUL value where the FIRST phase transition occurs
first_transition_points <- phase_transition_data %>%
  filter(is_phase_transition_step == 1) %>%
  group_by(primary_seed_node, social_distance_h) %>%
  summarise(
    first_transition_IUL = min(innovation_iul_Gamma),
    .groups = 'drop'
  )

# Now, for the heatmap: for each (IUL_sweep_value, h_sweep_value),
# count how many seeds had their first transition AT OR BEFORE that IUL_sweep_value.
heatmap_data_list <- list()
for (iul_val in IUL_values_sweep) {
  for (h_val in h_values_sweep) {
    count_transitions_at_or_before_IUL <- first_transition_points %>%
      filter(social_distance_h == h_val, first_transition_IUL <= iul_val) %>%
      nrow()
    
    prop_transitions <- count_transitions_at_or_before_IUL / length(primary_seed_nodes_for_study)
    
    heatmap_data_list[[length(heatmap_data_list) + 1]] <- data.frame(
      innovation_iul_Gamma = iul_val,
      social_distance_h = h_val,
      proportion_with_transition = prop_transitions
    )
  }
}
heatmap_df <- bind_rows(heatmap_data_list)

# -----------------------------------------------------------------------------
# 5. Visualization of Phase Transitions
# -----------------------------------------------------------------------------
cat("Generating phase transition heatmap...\n")

# Ensure h is treated as a factor for discrete y-axis, and Gamma as continuous
heatmap_df$social_distance_h_factor <- factor(sprintf("%.2f", heatmap_df$social_distance_h))

phase_transition_plot <- ggplot(heatmap_df, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = proportion_with_transition)) +
  geom_tile(color = "white", lwd = 0.2) + # Add white lines between tiles
  scale_fill_viridis_c(name = "Prop. of Seeds\nwith Phase Transition\n(Jump >= 0.5)", limits = c(0, 1), option="plasma") +
  labs(
    x = expression(paste("Intrinsic Innovation Utility (", Gamma, ")")),
    y = "Maximum Social Distance (h)",
    title = paste("Phase Transition Map for", CURRENT_GRAPH_TYPE_LABEL),
    subtitle = paste("Thresholds ~ N(μ=", THRESHOLD_MEAN_FIXED, ", σ=", TAU_NORMAL_DISTRIBUTION_SD_FIXED, "), ", length(primary_seed_nodes_for_study), " random seeds per (Γ,h) point", sep="")
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0,0)) +
  #scale_y_discrete(expand = c(0,0)) + # Keep if h values are not too many
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size=9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(phase_transition_plot)

# Save the plot
plot_filename <- paste0("phase_transition_heatmap_", gsub("[: ]", "_", Sys.time()), ".png")
ggsave(plot_filename, phase_transition_plot, width = 8, height = 6, dpi = 300)
cat(paste("Phase transition heatmap saved to:", plot_filename, "\n"))

# Save the raw data and heatmap data
saveRDS(all_sim_results_raw_df, "phase_transition_raw_results.rds")
saveRDS(heatmap_df, "phase_transition_heatmap_data.rds")
cat("Raw simulation results and heatmap data saved.\n")