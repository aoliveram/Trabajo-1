# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script: 06_GSS_diffusion_sims.R
#
# Objetivo:
#   Simula la difusión de un comportamiento (análogo a la adopción de una innovación)
#   sobre las 100 redes GSS previamente generadas.
#   El modelo de difusión es un modelo de umbral complejo que considera la
#   influencia social homofílica y la propensión individual a la acción.
#
# Genera (rds):
#   Archivos de resultados crudos en la carpeta de simulación, con nombres como:
#   trabajo_1_files/diffusion_simulation_files_sigm/GSS_phase_transition_raw_... .rds
#
# Fuentes:
#   - Redes 'GSS': trabajo_1_files/GSS_network_ergm/GSS_network_simulated_1000_mur_XXX.rds
#   - Funciones de simulación: diffusion_tests/simulation_functions.R
#
# Parámetros Clave:
#   - Los parámetros de simulación (umbrales, IUL, h) se mantienen idénticos
#     al script de ATP para permitir una comparación directa de los resultados
#     si se deseara en el futuro.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- 0. Cargar Librerías y Funciones ---
library(igraph)
library(doParallel)
library(dplyr)
library(network) # Necesario para asIgraph
library(intergraph)
library(cluster)

# Cargar el script que contiene las funciones de simulación
# Asegúrate de que la ruta sea correcta desde tu directorio de trabajo actual
source("diffusion_tests/simulation_functions.R")

# -----------------------------------------------------------------------------
# 1. Parámetros Principales y Configuración
# -----------------------------------------------------------------------------
cat("Configurando los parámetros principales de la simulación...\n")

# ---- MODIFICACIÓN 1: Cambiar etiquetas y directorios de ATP a GSS ----
CURRENT_GRAPH_TYPE_LABEL <- "GSS-net"
NETWORKS_DIR <- "trabajo_1_files/GSS_network_ergm/" # Directorio donde están las redes GSS con el score
# --------------------------------------------------------------------

# Parámetros (se mantienen sin cambios según lo solicitado)
NUM_NETWORK_INSTANCES_AVAILABLE <- 100 # Tenemos 100 redes GSS
NUM_SEED_RUNS_TOTAL <- 96 # Se usarán las primeras 96 redes para mantener consistencia con la corrida ATP
N_NODES_GLOBAL <- 1000

# Barrido de parámetros de umbral
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6)
TAU_NORMAL_SD_SWEEP_LIST <- c(0.08, 0.12, 0.16, 0.20)

# Configuración de paralelización
NUM_CORES_TO_USE <- 8 # M4

# Directorio de salida
RESULTS_DIR <- "trabajo_1_files/GSS_diffusion_simulation_files_sigm/"
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# 2. Parámetros de Barrido Interno (IUL y h)
# -----------------------------------------------------------------------------
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP <- seq(0/12, 12/12, by = 1/12)

# -----------------------------------------------------------------------------
# 3. Bucle Principal de Simulación
# -----------------------------------------------------------------------------
cat("Iniciando el gran barrido de simulación para redes GSS...\n")

# Seleccionar la estrategia de siembra. Cambia el valor aquí si es necesario.
strategies <- c("random", "central", "marginal", "eigen", "closeness")
SEEDING_STRATEGY_FIXED <- strategies[1] # Ejemplo: empezar con 'random'

all_sds_raw_results_list <- list()
total_sd_iterations <- length(TAU_NORMAL_SD_SWEEP_LIST)
time_init <- Sys.time()

for (sd_idx in 1:total_sd_iterations) {
  current_tau_sd <- TAU_NORMAL_SD_SWEEP_LIST[sd_idx]
  
  cat(paste0("\n####################################################################\n"))
  cat(paste0("Procesando para Desviación Estándar del Umbral (σ): ", current_tau_sd, "\n"))
  cat(paste0("####################################################################\n"))
  
  current_sd_all_means_raw_results <- list()
  total_mean_iterations <- length(THRESHOLD_MEAN_SWEEP_LIST)
  
  for (mean_idx in 1:total_mean_iterations) {
    current_threshold_mean <- THRESHOLD_MEAN_SWEEP_LIST[mean_idx]
    
    cat(paste0("\n====================================================================\n"))
    cat(paste0("  Procesando para Distribución del Umbral: Normal(μ=", current_threshold_mean, ", σ=", current_tau_sd, ")\n"))
    cat(paste0("====================================================================\n"))
    
    cl <- makeCluster(NUM_CORES_TO_USE, type = "FORK")
    registerDoParallel(cl)
    cat(paste0("    Cluster registrado con ", NUM_CORES_TO_USE, " cores.\n"))
    
    list_of_results_for_this_mean_sd_combo <- foreach(
      run_idx = 1:NUM_SEED_RUNS_TOTAL,
      .combine = 'list',
      .multicombine = TRUE,
      .packages = c('igraph', 'dplyr', 'intergraph', 'cluster'),
      .export = c('sweep_homoph_parameter', 'NETWORKS_DIR', 'current_threshold_mean', 
                  'current_tau_sd', 'IUL_VALUES_SWEEP', 'H_VALUES_SWEEP', 
                  'SEEDING_STRATEGY_FIXED', 'NUM_NETWORK_INSTANCES_AVAILABLE', 'N_NODES_GLOBAL'), 
      .errorhandling = 'pass'
    ) %dopar% {
      
      # ---- MODIFICACIÓN 2: Cambiar el nombre del archivo de red de ATP a GSS ----
      network_file_idx <- ((run_idx - 1) %% NUM_NETWORK_INSTANCES_AVAILABLE) + 1
      current_network_path <- file.path(NETWORKS_DIR, paste0("GSS_network_simulated_1000_mur_", sprintf("%03d", network_file_idx), ".rds"))
      # -----------------------------------------------------------------------
      
      if (!file.exists(current_network_path)) return(NULL) 
      
      graph_for_this_run_ergm <- readRDS(current_network_path)
      graph_for_this_run <- asIgraph(graph_for_this_run_ergm)
      N_NODES_SPECIFIC_GRAPH <- vcount(graph_for_this_run)
      
      # ---- MODIFICACIÓN 3: Usar 'propensity_score' en lugar de 'q_i' ----
      set.seed(run_idx * 3000 + round(current_threshold_mean * 100) + round(current_tau_sd * 100))
      node_mur_q_specific <- V(graph_for_this_run)$propensity_score
      # -------------------------------------------------------------------
      
      node_degrees_specific <- igraph::degree(graph_for_this_run)
      
      attributes_for_distance_specific <- data.frame(
        age = V(graph_for_this_run)$age, educ_num = V(graph_for_this_run)$educ_num,
        race = as.factor(V(graph_for_this_run)$race), relig = as.factor(V(graph_for_this_run)$relig),
        sex = as.factor(V(graph_for_this_run)$sex)
      )
      d_ij_matrix <- as.matrix(daisy(attributes_for_distance_specific, metric = "gower"))
      
      set.seed(run_idx * 1000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000)) 
      node_thresholds_tau_frac_specific <- rnorm(n = N_NODES_SPECIFIC_GRAPH, mean = current_threshold_mean, sd = current_tau_sd)
      node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific < 0] <- 0
      node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific > 1] <- 1
      
      node_thresholds_count_for_cluster_specific <- round(node_thresholds_tau_frac_specific * node_degrees_specific)
      node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific > 0] <- 1
      node_thresholds_count_for_cluster_specific[node_thresholds_count_for_cluster_specific <= 0 & node_thresholds_tau_frac_specific == 0] <- 0
      
      # Lógica de siembra (sin cambios, usa la estrategia definida)
      if (SEEDING_STRATEGY_FIXED == "central") {
        primary_seed_for_this_run <- which.max(igraph::degree(graph_for_this_run))
      } else if (SEEDING_STRATEGY_FIXED == "marginal") {
        lowest_10_percent <- sort(igraph::degree(graph_for_this_run), index.return = TRUE)$ix[1:ceiling(N_NODES_GLOBAL * 0.1)]
        set.seed(run_idx * 2000)
        primary_seed_for_this_run <- as.integer(sample(lowest_10_percent, 1))
      } else if (SEEDING_STRATEGY_FIXED == "closeness") {
        primary_seed_for_this_run <- which.max(igraph::closeness(graph_for_this_run))
      } else if (SEEDING_STRATEGY_FIXED == "eigen") {
        primary_seed_for_this_run <- which.max(igraph::eigen_centrality(graph_for_this_run)$vector)
      } else { # "random" es el default
        set.seed(run_idx * 2000)
        primary_seed_for_this_run <- as.integer(sample(V(graph_for_this_run), 1))
      }
      if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
      
      num_seeds_for_initial_cluster <- node_thresholds_count_for_cluster_specific[primary_seed_for_this_run]
      num_seeds_for_initial_cluster <- min(num_seeds_for_initial_cluster, N_NODES_SPECIFIC_GRAPH, node_degrees_specific[primary_seed_for_this_run] + 1, na.rm = TRUE)
      if (num_seeds_for_initial_cluster < 1) num_seeds_for_initial_cluster <- 1
      
      initial_infectors_for_this_sim_run <- c(primary_seed_for_this_run)
      if (num_seeds_for_initial_cluster > 1) {
        neighbors_of_primary <- as.numeric(neighbors(graph_for_this_run, primary_seed_for_this_run, mode="all"))
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
      
      # Llamada a la función de simulación (sin cambios en sus argumentos)
      df_one_full_run <- sweep_homoph_parameter(
        primary_seed_id_arg = primary_seed_for_this_run, 
        N_nodes_arg = N_NODES_SPECIFIC_GRAPH, 
        graph_obj_arg = graph_for_this_run,
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
    } # Fin del loop foreach
    
    stopCluster(cl)
    cat(paste0("    Simulaciones paralelas para (μ=", current_threshold_mean, ", σ=", current_tau_sd, ") finalizadas.\n"))
    
    valid_results <- !sapply(list_of_results_for_this_mean_sd_combo, function(x) inherits(x, "simpleError") || is.null(x))
    if(sum(valid_results) > 0) {
      current_mean_sd_raw_df <- bind_rows(list_of_results_for_this_mean_sd_combo[valid_results])
      current_sd_all_means_raw_results[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- current_mean_sd_raw_df
    } else {
      cat(paste0("    ADVERTENCIA: No se obtuvieron resultados válidos para esta combinación de parámetros.\n"))
      current_sd_all_means_raw_results[[paste0("mean_", sprintf("%.2f", current_threshold_mean))]] <- NULL
    }
    
  } # Fin del loop de 'mean'
  
  # ---- MODIFICACIÓN 4: Cambiar el nombre del archivo de salida de ATP a GSS ----
  raw_data_filename <- file.path(RESULTS_DIR, paste0("GSS_phase_transition_raw_sd", sprintf("%.2f", current_tau_sd), "_means_all_", SEEDING_STRATEGY_FIXED, ".rds"))
  saveRDS(current_sd_all_means_raw_results, raw_data_filename)
  cat(paste0("    Resultados crudos para σ=", current_tau_sd, " guardados en: ", raw_data_filename, "\n"))
  # -------------------------------------------------------------------------
  
  all_sds_raw_results_list[[paste0("sd_", sprintf("%.2f", current_tau_sd))]] <- current_sd_all_means_raw_results
  
} # Fin del loop de 'sd'

time_fin <- Sys.time()
time_total <- difftime(time_fin, time_init, units = "auto")
cat(paste("\nTodos los barridos de simulación han finalizado. Tiempo total:", format(time_total), "\n"))

# -----------------------------------------------------------------------------
# 4. Guardado Final de Datos Agregados
# -----------------------------------------------------------------------------
cat("Guardando el objeto combinado final de resultados...\n")

# ---- MODIFICACIÓN 5: Cambiar el nombre del archivo de salida final de ATP a GSS ----
final_output_filename <- file.path(RESULTS_DIR, paste0("GSS_phase_transition_GRAND_COMBINED_raw_results_", SEEDING_STRATEGY_FIXED, ".rds"))
saveRDS(all_sds_raw_results_list, final_output_filename)
# ----------------------------------------------------------------------------

cat(paste("Todos los datos han sido guardados. El script ha finalizado. Archivo principal:", final_output_filename, "\n"))