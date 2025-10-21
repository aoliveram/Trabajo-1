# =============================================================================
# 04_GSS_diffusion_sims_Slurm.R
# Run with Slurm (PSOCK multinode). Example sbatch:
#
# #!/bin/bash
# #SBATCH --job-name=gss-psock
# #SBATCH --account=vegayon-np
# #SBATCH --partition=vegayon-np
# #SBATCH --nodes=4
# #SBATCH --ntasks-per-node=1
# #SBATCH --cpus-per-task=16
# #SBATCH --time=08:00:00
# #SBATCH --output=slurm-%j.out
# #SBATCH --error=slurm-%j.err
# export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# module load R/4.4.2
# srun Rscript 04_GSS_diffusion_sims_Slurm.R
# =============================================================================

# --- 0) Libraries & helpers ---------------------------------------------------
suppressPackageStartupMessages({
  library(igraph)
  library(doParallel)
  library(dplyr)
  library(network)     # for asIgraph (via intergraph)
  library(intergraph)
  library(cluster)     # daisy (Gower)
  library(foreach)
})

# Reproducibility knobs for BLAS/OpenMP per worker
Sys.setenv(OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1")

# Source simulation functions (master); workers will source this too
source("diffusion_tests/simulation_functions.R")

# --- 1) Parameters & configuration -------------------------------------------
cat("Configurando parámetros principales...\n")

# **Always use ER** as requested
CURRENT_GRAPH_TYPE_LABEL <- "ER"

# Paths
NETWORKS_DIR <- "trabajo_1_files/GSS_network_ergm/"
RESULTS_DIR  <- "trabajo_1_files/GSS_ER_diffusion_simulation_files_sigm/"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# Data availability / size
NUM_NETWORK_INSTANCES_AVAILABLE <- 100
NUM_SEED_RUNS_TOTAL             <- 96
N_NODES_GLOBAL                  <- 1000

# Sweeps (outer: means, sds)
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6)
TAU_NORMAL_SD_SWEEP_LIST  <- c(0.08, 0.12, 0.16, 0.20)

# Inner sweeps (IUL, h)
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP   <- seq(0/12, 12/12, by = 1/12)

# Seed strategies to run (you already ran "random")
strategies <- c("central", "marginal", "eigen", "closeness")

# --- 2) PSOCK multinode cluster from Slurm allocation ------------------------
# Detect Slurm hosts and build worker list
get_slurm_hosts <- function() {
  h <- tryCatch(system("scontrol show hostnames $SLURM_JOB_NODELIST", intern = TRUE),
                error = function(e) character(0))
  if (length(h) == 0) "localhost" else h
}
SLURM_HOSTS    <- get_slurm_hosts()
CORES_PER_NODE <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
worker_hosts   <- rep(SLURM_HOSTS, each = CORES_PER_NODE)

cat("Slurm hosts detectados:\n")
print(SLURM_HOSTS)
cat("Cores por nodo:", CORES_PER_NODE, "\n")
cat("Workers totales:", length(worker_hosts), "\n")

# Create PSOCK cluster once and reuse for all loops
cl <- parallel::makeCluster(worker_hosts, type = "PSOCK", outfile = "")
doParallel::registerDoParallel(cl)

# Ensure required packages and functions are present on workers
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(igraph); library(dplyr); library(intergraph); library(cluster)
  })
  # Pin thread envs in workers as well
  Sys.setenv(OMP_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1",
             OPENBLAS_NUM_THREADS = "1")
  NULL
})

# Workers need helpers (either export functions or source file on each worker)
# Sourcing is simpler/robust because functions live there:
parallel::clusterExport(cl, c("NETWORKS_DIR", "IUL_VALUES_SWEEP", "H_VALUES_SWEEP",
                              "NUM_NETWORK_INSTANCES_AVAILABLE", "N_NODES_GLOBAL"),
                        envir = environment())
parallel::clusterEvalQ(cl, {
  source("diffusion_tests/simulation_functions.R")
  NULL
})

# Helpers defined in THIS file (workers need them too)
graph_density_target <- function(g) {
  igraph::ecount(g) / (igraph::vcount(g) * (igraph::vcount(g) - 1) / 2)
}
random_er_same_density <- function(g_base) {
  n <- igraph::vcount(g_base)
  p <- graph_density_target(g_base)
  gr <- igraph::erdos.renyi.game(n = n, p.or.m = p, type = "gnp", directed = FALSE, loops = FALSE)
  if (!is.null(igraph::V(g_base)$name)) igraph::V(gr)$name <- igraph::V(g_base)$name
  for (attr in igraph::vertex_attr_names(g_base)) {
    igraph::vertex_attr(gr, attr) <- igraph::vertex_attr(g_base, attr)
  }
  gr
}
random_degree_preserving <- function(g_base, niter_factor = 20) {
  m  <- igraph::ecount(g_base)
  gr <- igraph::rewire(g_base, with = igraph::keeping_degseq(niter = niter_factor * m))
  gr
}
parallel::clusterExport(cl, c("graph_density_target",
                              "random_er_same_density",
                              "random_degree_preserving"),
                        envir = environment())

# --- 3) Main loops -----------------------------------------------------------
cat("Iniciando barrido de simulación GSS (topología ER) ...\n")
time_init <- Sys.time()

for (SEEDING_STRATEGY_FIXED in strategies) {
  cat("\n================ STRATEGY:", SEEDING_STRATEGY_FIXED, "================\n")
  
  all_sds_raw_results_list <- list()
  total_sd_iterations      <- length(TAU_NORMAL_SD_SWEEP_LIST)
  
  for (sd_idx in seq_len(total_sd_iterations)) {
    current_tau_sd <- TAU_NORMAL_SD_SWEEP_LIST[sd_idx]
    
    cat("\n####################################################################\n")
    cat("σ (SD del umbral) =", current_tau_sd, "\n")
    cat("####################################################################\n")
    
    current_sd_all_means_raw_results <- list()
    total_mean_iterations            <- length(THRESHOLD_MEAN_SWEEP_LIST)
    
    for (mean_idx in seq_len(total_mean_iterations)) {
      current_threshold_mean <- THRESHOLD_MEAN_SWEEP_LIST[mean_idx]
      
      cat("\n====================================================================\n")
      cat("Normal(μ =", current_threshold_mean, ", σ =", current_tau_sd, ")\n")
      cat("====================================================================\n")
      
      # PSOCK cluster already registered; run foreach in parallel across ALL workers
      list_of_results_for_this_mean_sd_combo <- foreach::foreach(
        run_idx = 1:NUM_SEED_RUNS_TOTAL,
        .combine = 'list',
        .multicombine = TRUE,
        .maxcombine = NUM_SEED_RUNS_TOTAL,
        .packages = c('igraph', 'dplyr', 'intergraph', 'cluster'),
        .export   = c('sweep_homoph_parameter', 'NETWORKS_DIR', 'current_threshold_mean',
                      'current_tau_sd', 'IUL_VALUES_SWEEP', 'H_VALUES_SWEEP',
                      'SEEDING_STRATEGY_FIXED', 'NUM_NETWORK_INSTANCES_AVAILABLE',
                      'N_NODES_GLOBAL', 'CURRENT_GRAPH_TYPE_LABEL',
                      'graph_density_target', 'random_er_same_density',
                      'random_degree_preserving'),
        .errorhandling = 'pass'
      ) %dopar% {
        
        # Select network file
        network_file_idx     <- ((run_idx - 1) %% NUM_NETWORK_INSTANCES_AVAILABLE) + 1
        current_network_path <- file.path(NETWORKS_DIR,
                                          paste0("GSS_network_simulated_1000_mur_",
                                                 sprintf("%03d", network_file_idx), ".rds"))
        if (!file.exists(current_network_path)) return(NULL)
        
        graph_for_this_run_ergm <- readRDS(current_network_path)
        graph_for_this_run      <- asIgraph(graph_for_this_run_ergm)
        
        # --- Topology selection (always ER here) ---
        base_graph_for_attributes <- graph_for_this_run
        if (CURRENT_GRAPH_TYPE_LABEL == "ER") {
          graph_for_this_run <- random_er_same_density(base_graph_for_attributes)
        } else if (CURRENT_GRAPH_TYPE_LABEL == "ER_degseq") {
          graph_for_this_run <- random_degree_preserving(base_graph_for_attributes, niter_factor = 20)
          if (!is.null(igraph::V(base_graph_for_attributes)$name)) {
            igraph::V(graph_for_this_run)$name <- igraph::V(base_graph_for_attributes)$name
          }
          for (attr in igraph::vertex_attr_names(base_graph_for_attributes)) {
            if (!(attr %in% igraph::vertex_attr_names(graph_for_this_run))) {
              igraph::vertex_attr(graph_for_this_run, attr) <- igraph::vertex_attr(base_graph_for_attributes, attr)
            }
          }
        } else {
          # Keep original GSS topology
        }
        
        N_NODES_SPECIFIC_GRAPH <- igraph::vcount(graph_for_this_run)
        
        # Node attributes and distances
        set.seed(run_idx * 3000 + round(current_threshold_mean * 100) + round(current_tau_sd * 100))
        node_mur_q_specific <- igraph::V(graph_for_this_run)$propensity_score
        
        node_degrees_specific <- igraph::degree(graph_for_this_run)
        
        attributes_for_distance_specific <- data.frame(
          age      = igraph::V(graph_for_this_run)$age,
          educ_num = igraph::V(graph_for_this_run)$educ_num,
          race     = as.factor(igraph::V(graph_for_this_run)$race),
          relig    = as.factor(igraph::V(graph_for_this_run)$relig),
          sex      = as.factor(igraph::V(graph_for_this_run)$sex)
        )
        d_ij_matrix <- as.matrix(cluster::daisy(attributes_for_distance_specific, metric = "gower"))
        
        # Thresholds τ_i
        set.seed(run_idx * 1000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000))
        node_thresholds_tau_frac_specific <- rnorm(
          n    = N_NODES_SPECIFIC_GRAPH,
          mean = current_threshold_mean,
          sd   = current_tau_sd
        )
        node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific < 0] <- 0
        node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific > 1] <- 1
        
        node_thresholds_count_for_cluster_specific <- round(node_thresholds_tau_frac_specific * node_degrees_specific)
        node_thresholds_count_for_cluster_specific[
          node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific > 0] <- 1
        node_thresholds_count_for_cluster_specific[
          node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific == 0] <- 0
        
        # Seeding logic
        if (SEEDING_STRATEGY_FIXED == "central") {
          primary_seed_for_this_run <- which.max(igraph::degree(graph_for_this_run))
        } else if (SEEDING_STRATEGY_FIXED == "marginal") {
          lowest_10_percent <- sort(igraph::degree(graph_for_this_run),
                                    index.return = TRUE)$ix[1:ceiling(N_NODES_GLOBAL * 0.1)]
          set.seed(run_idx * 2000)
          primary_seed_for_this_run <- as.integer(sample(lowest_10_percent, 1))
        } else if (SEEDING_STRATEGY_FIXED == "closeness") {
          primary_seed_for_this_run <- which.max(igraph::closeness(graph_for_this_run))
        } else if (SEEDING_STRATEGY_FIXED == "eigen") {
          primary_seed_for_this_run <- which.max(igraph::eigen_centrality(graph_for_this_run)$vector)
        } else {
          set.seed(run_idx * 2000)
          primary_seed_for_this_run <- as.integer(sample(igraph::V(graph_for_this_run), 1))
        }
        if (length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
        
        num_seeds_for_initial_cluster <- node_thresholds_count_for_cluster_specific[primary_seed_for_this_run]
        num_seeds_for_initial_cluster <- min(num_seeds_for_initial_cluster,
                                             N_NODES_SPECIFIC_GRAPH,
                                             node_degrees_specific[primary_seed_for_this_run] + 1,
                                             na.rm = TRUE)
        if (num_seeds_for_initial_cluster < 1) num_seeds_for_initial_cluster <- 1
        
        initial_infectors_for_this_sim_run <- c(primary_seed_for_this_run)
        if (num_seeds_for_initial_cluster > 1) {
          neighbors_of_primary  <- as.numeric(igraph::neighbors(graph_for_this_run, primary_seed_for_this_run, mode = "all"))
          num_additional_needed <- num_seeds_for_initial_cluster - 1
          if (length(neighbors_of_primary) >= num_additional_needed) {
            initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run,
                                                    sample(neighbors_of_primary, num_additional_needed))
          } else {
            initial_infectors_for_this_sim_run <- c(initial_infectors_for_this_sim_run, neighbors_of_primary)
            still_needed_more <- num_seeds_for_initial_cluster - length(initial_infectors_for_this_sim_run)
            if (still_needed_more > 0) {
              potential_others <- setdiff(1:N_NODES_SPECIFIC_GRAPH, initial_infectors_for_this_sim_run)
              if (length(potential_others) > 0) {
                initial_infectors_for_this_sim_run <- c(
                  initial_infectors_for_this_sim_run,
                  sample(potential_others, min(length(potential_others), still_needed_more))
                )
              }
            }
          }
        }
        initial_infectors_for_this_sim_run <- unique(initial_infectors_for_this_sim_run)
        
        # Run simulation sweep for this run
        df_one_full_run <- sweep_homoph_parameter(
          primary_seed_id_arg              = primary_seed_for_this_run,
          N_nodes_arg                      = N_NODES_SPECIFIC_GRAPH,
          graph_obj_arg                    = graph_for_this_run,
          node_individual_thresholds_tau_arg = node_thresholds_tau_frac_specific,
          node_mur_q_arg                   = node_mur_q_specific,
          all_innovation_iul_Gamma_values  = IUL_VALUES_SWEEP,
          all_social_distance_h_values     = H_VALUES_SWEEP,
          initial_infectors_vector_arg     = initial_infectors_for_this_sim_run,
          d_ij_matrix                      = d_ij_matrix
        )
        
        df_one_full_run$run_id                    <- run_idx
        df_one_full_run$network_instance_file_idx <- network_file_idx
        df_one_full_run$threshold_mean_param      <- current_threshold_mean
        df_one_full_run$threshold_sd_param        <- current_tau_sd
        df_one_full_run$N_nodes_actual            <- N_NODES_SPECIFIC_GRAPH
        
        df_one_full_run
      } # %dopar%
      
      # Collect valid results
      valid_results <- !sapply(list_of_results_for_this_mean_sd_combo,
                               function(x) inherits(x, "simpleError") || is.null(x))
      if (sum(valid_results) > 0) {
        current_mean_sd_raw_df <- bind_rows(list_of_results_for_this_mean_sd_combo[valid_results])
        current_sd_all_means_raw_results[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- current_mean_sd_raw_df
      } else {
        cat("    ADVERTENCIA: sin resultados válidos para esta combinación μ,σ.\n")
        current_sd_all_means_raw_results[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- NULL
      }
    } # means loop
    
    # Save per-σ raw results (include strategy)
    raw_data_filename <- file.path(
      RESULTS_DIR,
      paste0("GSS_phase_transition_raw_sd", sprintf("%.2f", current_tau_sd),
             "_means_all_", SEEDING_STRATEGY_FIXED, ".rds")
    )
    saveRDS(current_sd_all_means_raw_results, raw_data_filename)
    cat("    Guardado resultados σ =", current_tau_sd, "→", raw_data_filename, "\n")
    
    all_sds_raw_results_list[[paste0("sd_", sprintf("%.2f", current_tau_sd))]] <- current_sd_all_means_raw_results
  } # sds loop
  
  # Save combined per-strategy
  final_output_filename <- file.path(
    RESULTS_DIR,
    paste0("GSS_phase_transition_GRAND_COMBINED_raw_results_", SEEDING_STRATEGY_FIXED, ".rds")
  )
  saveRDS(all_sds_raw_results_list, final_output_filename)
  cat("Archivo combinado por strategy guardado →", final_output_filename, "\n")
  
} # strategies loop

time_fin   <- Sys.time()
time_total <- difftime(time_fin, time_init, units = "auto")
cat("\nTodos los barridos han finalizado. Tiempo total:", format(time_total), "\n")

# --- 4) Clean up --------------------------------------------------------------
stopCluster(cl)
cat("Cluster PSOCK detenido. Fin del script.\n")