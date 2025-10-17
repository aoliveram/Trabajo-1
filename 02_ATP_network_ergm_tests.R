# ===============================================================
# BATERÍA DE TESTS TOPOLOGÍA REDES ATP (N=1000), 100 realizaciones
# ---------------------------------------------------------------
# Entradas:
#   - Directorio con RDS de redes simuladas (clase 'network' de {ergm})
# Salidas (todas comentadas; descomentar para escribir):
#   - results/*.csv    trabajo_1_plots/ATP_network_tests/*.png
# ===============================================================


library(igraph)
library(Matrix)
library(dplyr)
library(purrr)
library(readr)
library(tibble)
library(stringr)

library(ggplot2)
library(network)  # needed for as.matrix.network()
library(tidyr)    # needed for pivot_longer()

library(parallel)
library(doParallel)

library(foreach)
library(patchwork)    # para combinar plots (lineal + log-log en una sola figura)

dir_input <- "trabajo_1_files/ATP_network_ergm"
pattern_rds <- "^ATP_network_simulated_1000_\\d{3}\\.rds$"

# ==============================================================================
# Estadísticos
# ==============================================================================


# Control de tamaños de submuestreo para small-world (Estrategia A)
subsample_sizes <- c(200, 400, 600, 800, 1000)
subsample_reps  <- 5  # repeticiones por tamaño

# Rewiring bajo (fracción de aristas a re cablear)
rewire_frac <- 0.01

# Pasos de percolación (fracción de nodos removidos)
percolation_p <- seq(0, 1, by = 0.02)

# Semilla reproducible para submuestreos/rewire
set.seed(12345)

# ---------- Utilidades ----------
# Cargar 'network' (ergm) -> igraph simple no dirigido
network_rds_to_igraph <- function(path) {
  net <- readRDS(path)
  # Extraer matriz de adyacencia desde 'network'
  adj <- as.matrix.network(net, matrix.type = "adjacency")
  g   <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  g
}

giant_component <- function(g) {
  comps <- components(g)
  induced_subgraph(g, vids = which(comps$membership == which.max(comps$csize)))
}

avg_path_and_diameter <- function(g) {
  if (ecount(g) == 0 || vcount(g) < 2) return(c(L = NA_real_, D = NA_real_))
  G <- giant_component(g)
  c(L = suppressWarnings(average.path.length(G, unconnected = FALSE)),
    D = suppressWarnings(diameter(G, directed = FALSE, unconnected = FALSE)))
}

assort_deg <- function(g) {
  if (ecount(g) == 0) return(NA_real_)
  assortativity_degree(g, directed = FALSE)
}

clustering_all <- function(g) {
  c(C_global = transitivity(g, type = "global"),
    C_avg    = mean(transitivity(g, type = "local", isolates = "zero"), na.rm = TRUE))
}

gini_vec <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  G <- (2 * sum(x * seq_len(n)))/(n * sum(x)) - (n + 1)/n
  G
}

# Índice small-world de Humphries & Gurney (2008): sigma = (C/C_rand) / (L/L_rand)
smallworld_sigma <- function(g, g_rand) {
  L <- avg_path_and_diameter(g)["L"]
  C <- transitivity(g, "global")
  Lr <- avg_path_and_diameter(g_rand)["L"]
  Cr <- transitivity(g_rand, "global")
  if (any(is.na(c(L, C, Lr, Cr))) || Lr == 0 || Cr == 0) return(NA_real_)
  (C / Cr) / (L / Lr)
}

# ER con N y m exactos
er_nm <- function(N, m) {
  g <- sample_gnm(N, m, directed = FALSE, loops = FALSE)
  igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

# Watts–Strogatz con mismo N y grado medio aprox. (k par)
ws_match <- function(N, m, p = 0.05) {
  k <- round(2 * m / N)
  if (k %% 2 == 1) k <- k + 1
  if (k < 2) k <- 2
  g <- igraph::sample_smallworld(dim = 1, size = N, nei = k/2, p = p, loops = FALSE, multiple = FALSE)
  igraph::simplify(g)
}

# Barabási–Albert (BA): escoger 'm0' tal que <k> ~ 2m0 ≈ 2m/N
ba_match <- function(N, m) {
  mm <- max(1, round(m / N))  # edges added per step
  g <- sample_pa(n = N, power = 1, m = mm, directed = FALSE)
  igraph::simplify(g)
}

# Configuration model (misma secuencia de grados); simplificar multiaristas/bucles
config_match <- function(deg) {
  deg <- as.integer(deg)
  # Asegurar que la secuencia sea "graphical" y factible
  if (sum(deg) %% 2 == 1) {
    # Forzar suma par ajustando el nodo con mayor grado
    idx <- which.max(deg)
    deg[idx] <- deg[idx] + 1L
  }
  # Evitar grados inválidos (>= N)
  N <- length(deg)
  deg[deg >= N] <- N - 1L
  # Intento 1: construcción simple sin multiaristas/loops (Havel-Hakimi)
  g <- tryCatch(
    igraph::sample_degseq(deg, method = "simple.no.multiple"),
    error = function(e) NULL
  )
  # Intento 2: construcción simple permitiendo multiaristas, luego simplificar
  if (is.null(g)) {
    g <- tryCatch(
      igraph::sample_degseq(deg, method = "simple"),
      error = function(e) NULL
    )
    if (!is.null(g)) {
      g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    }
  }
  # Intento 3: último recurso — configuration model por emparejamiento ("vl"); puede fallar
  if (is.null(g)) {
    g <- tryCatch(
      igraph::sample_degseq(deg, method = "vl"),
      error = function(e) NULL
    )
    if (!is.null(g)) {
      g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    }
  }
  if (is.null(g)) {
    stop("No se pudo construir un grafo simple para la secuencia de grados dada.")
  }
  g
}

# Rewire fracción de aristas: eliminar p*m aristas al azar y reconectar extremos a no-enlaces
rewire_fraction <- function(g, frac = 0.01) {
  if (ecount(g) == 0 || frac <= 0) return(g)
  # Maslov–Sneppen edge swaps that preserve the degree sequence
  niter <- max(1L, as.integer(round(frac * ecount(g)) * 2L))
  igraph::rewire(g, with = keeping_degseq(niter = niter))
}

# Small-world por submuestreo (Estrategia A)
L_vs_logN_subsampling <- function(g, sizes = subsample_sizes, reps = subsample_reps) {
  N <- vcount(g)
  sizes <- sort(unique(pmin(sizes, N)))
  out <- list()
  vids_all <- V(g)
  for (n_sub in sizes) {
    Lvals <- numeric(reps)
    for (r in seq_len(reps)) {
      idx <- sample(vids_all, n_sub)
      sg  <- induced_subgraph(g, idx)
      Lvals[r] <- avg_path_and_diameter(sg)["L"]
    }
    out[[as.character(n_sub)]] <- data.frame(Nprime = n_sub,
                                             logN   = log(n_sub),
                                             L      = median(Lvals, na.rm = TRUE))
  }
  bind_rows(out)
}

# Triadic census en no dirigidas (0,1,2,3 aristas por triada)
# Cálculo eficiente usando matriz dispersa y triángulos
triadic_census_undirected <- function(g) {
  N  <- vcount(g)
  Tt <- choose(N, 3)  # total triadas
  deg <- degree(g)
  # Triángulos
  tri <- count_triangles(g)          # por-vertice
  n_tri <- sum(tri) / 3              # total triángulos (cada triángulo contado 3 veces)
  # Wedges (caminos de longitud 2): sum choose(deg,2) - 3*triángulos
  W <- sum(choose(deg, 2)) - 3 * n_tri
  # Para triadas de 1 arista: S = sum_{(u,v)∈E} [N - deg(u) - deg(v) + c_uv]
  # donde c_uv = # vecinos comunes = (A^2)_{uv}
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  A2 <- A %*% A
  Edf <- as.data.frame(as_edgelist(g, names = FALSE))
  colnames(Edf) <- c("u","v")
  # extraer c_uv desde A2 (nota: A2 simétrica)
  get_c_uv <- function(u, v) A2[u, v]
  c_uv <- mapply(get_c_uv, Edf$u, Edf$v)
  S <- sum(N - deg[Edf$u] - deg[Edf$v] + c_uv)
  # 0 aristas: Z = T - (S + W + n_tri)
  Z <- Tt - (S + W + n_tri)
  tibble(
    triads_0 = as.numeric(Z),
    triads_1 = as.numeric(S),
    triads_2 = as.numeric(W),
    triads_3 = as.numeric(n_tri)
  )
}

# C(k) vs k (promedio local por grado)
Ck_curve <- function(g) {
  deg <- degree(g)
  Cl  <- transitivity(g, type = "local", isolates = "zero")
  df  <- tibble(k = deg, C = Cl) %>%
    filter(is.finite(k), is.finite(C)) %>%
    group_by(k) %>% summarise(Ck = mean(C, na.rm = TRUE), n = n(), .groups="drop")
  df
}

# Percolación: tamaño GCC bajo fallas aleatorias / ataques dirigidos
percolation_curves <- function(g, p_grid = percolation_p, mode = c("random", "attack"), re_rank_static = TRUE) {
  mode <- match.arg(mode)
  N <- vcount(g)
  baseline_deg <- degree(g)
  order_attack <- order(baseline_deg, decreasing = TRUE)
  S <- numeric(length(p_grid))
  for (i in seq_along(p_grid)) {
    p <- p_grid[i]
    n_remove <- round(p * N)
    if (mode == "random") {
      rem <- if (n_remove > 0) sample.int(N, n_remove) else integer(0)
    } else {
      if (re_rank_static) {
        rem <- order_attack[seq_len(min(n_remove, N))]
      } else {
        # re-rank dinámico: más costoso; opcional
        tmpg <- g
        rem <- integer(0)
        for (k in seq_len(n_remove)) {
          d <- degree(tmpg)
          if (length(d) == 0) break
          v <- which.max(d)
          rem <- c(rem, as.integer(V(tmpg)[v]))
          tmpg <- delete_vertices(tmpg, v)
        }
      }
    }
    g_sub <- if (length(rem) > 0) delete_vertices(g, rem) else g
    if (vcount(g_sub) == 0 || ecount(g_sub) == 0) {
      S[i] <- 0
    } else {
      S[i] <- components(g_sub)$csize %>% max() %>% `/`(N) %>% as.numeric()
    }
  }
  tibble(p = p_grid, S = S, mode = mode)
}

# ---------- Carga de las 100 redes ----------
files <- list.files(dir_input, pattern = pattern_rds, full.names = TRUE)

# Backend paralelo con FORK
cl <- makeCluster(8, type = "FORK")
registerDoParallel(cl)

process_one_graph <- function(path_rds) {
  g <- network_rds_to_igraph(path_rds)
  N <- vcount(g); m <- ecount(g)
  L_D <- avg_path_and_diameter(g)
  C_all <- clustering_all(g)
  base <- c(
    N       = N,
    m       = m,
    density = edge_density(g),
    k_mean  = mean(degree(g)),
    L       = unname(L_D["L"]),
    D       = unname(L_D["D"]),
    assort_deg = assort_deg(g),
    C_global   = unname(C_all["C_global"]),
    C_avg      = unname(C_all["C_avg"]),
    gini_deg   = gini_vec(degree(g))
  )
  # Small-world por submuestreo
  sw_sub <- L_vs_logN_subsampling(g)
  # Índice small-world (sigma) vs ER empírico
  g_er <- er_nm(N, m)
  sigma <- smallworld_sigma(g, g_er)
  # Rewiring bajo
  g_rw <- rewire_fraction(g, frac = rewire_frac)
  L_rw_val <- avg_path_and_diameter(g_rw)["L"]
  if (is.na(L_rw_val)) {
    Gtmp <- giant_component(g_rw)
    L_rw_val <- if (vcount(Gtmp) >= 2 && ecount(Gtmp) > 0) mean_distance(Gtmp, directed = FALSE, unconnected = FALSE) else NA_real_
  }
  rw_metrics <- c(L_rw = as.numeric(L_rw_val),
                  C_rw = transitivity(g_rw, "global"),
                  rdeg_rw = assort_deg(g_rw))
  # Modelos de referencia para comparación de grados
  deg_g <- degree(g)
  g_ws <- ws_match(N, m, p = 0.05)
  g_ba <- ba_match(N, m)
  g_cf <- config_match(deg_g)
  deg_tbl <- tibble(
    model = c("ATP","ER","WS","BA","CONFIG"),
    gini  = c(gini_vec(deg_g),
              gini_vec(degree(g_er)),
              gini_vec(degree(g_ws)),
              gini_vec(degree(g_ba)),
              gini_vec(degree(g_cf))),
    cv    = c(sd(deg_g)/mean(deg_g),
              sd(degree(g_er))/mean(degree(g_er)),
              sd(degree(g_ws))/mean(degree(g_ws)),
              sd(degree(g_ba))/mean(degree(g_ba)),
              sd(degree(g_cf))/mean(degree(g_cf)))
  )
  # Censo triádico (no dirigido)
  tri_tbl <- triadic_census_undirected(g)
  # C(k) vs k
  ck_tbl  <- Ck_curve(g)
  # Percolación
  perc_random <- percolation_curves(g, mode = "random")
  perc_attack <- percolation_curves(g, mode = "attack")
  list(
    id   = basename(path_rds),
    g    = g,
    base = base,
    sw_sub = sw_sub,
    sigma = sigma,
    rewiring = rw_metrics,
    degree_cmp = deg_tbl,
    triad = tri_tbl,
    Ck = ck_tbl,
    percolation = bind_rows(perc_random, perc_attack) %>% mutate(id = basename(path_rds))
  )
}

message("Cargando y procesando ", length(files), " redes...")

res_list <- foreach(
  file = files,
  .export = c(
    "process_one_graph",
    "network_rds_to_igraph","giant_component","avg_path_and_diameter",
    "assort_deg","clustering_all","gini_vec","smallworld_sigma",
    "er_nm","ws_match","ba_match","config_match","rewire_fraction",
    "L_vs_logN_subsampling","triadic_census_undirected","Ck_curve","percolation_curves"
  ),
  .packages = c("igraph","Matrix","dplyr","purrr","readr","tibble","stringr","ggplot2","tidyr","network")
) %dopar% {
  process_one_graph(file)
}
stopCluster(cl)
registerDoSEQ()

# ---------- Ensambles / “condensados” ----------
# Tabla 1: métricas base por red
metrics_by_net <- map_dfr(res_list, ~{
  as_tibble_row(.x$base) %>% mutate(id = .x$id, .before = 1)
})

# Resumen agregado (media, mediana, sd, IQR)
metrics_summary <- metrics_by_net %>%
  select(-id) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  summarise(
    mean    = mean(value, na.rm = TRUE),
    median  = median(value, na.rm = TRUE),
    sd      = sd(value, na.rm = TRUE),
    p25     = quantile(value, 0.25, na.rm = TRUE),
    p75     = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Número de componentes por red
# La mitad son conexas. El resto tiene 1 o 2 aislados. 
components_by_net <- map_dfr(res_list, ~{
  comps <- igraph::components(.x$g)
  sizes <- sort(as.integer(comps$csize), decreasing = TRUE)
  tibble(
    id           = .x$id,
    n_components = length(sizes),
    gcc_size     = sizes[1],
    sizes        = list(sizes)
  )
})

# Tamaño medio de componentes pequeñas cuando n_components > 1
small_comps_by_net <- components_by_net %>%
  filter(n_components > 1) %>%
  transmute(
    id,
    small_sizes = purrr::map(sizes, ~ if (length(.x) >= 2) .x[-1] else integer(0)),
    n_small     = purrr::map_int(small_sizes, length),
    mean_small  = purrr::map_dbl(small_sizes, ~ if (length(.x) > 0) mean(.x) else NA_real_),
    med_small   = purrr::map_dbl(small_sizes, ~ if (length(.x) > 0) median(.x) else NA_real_),
    prop_isol   = purrr::map_dbl(small_sizes, ~ if (length(.x) > 0) mean(.x == 1) else NA_real_)
  )

small_comps_summary <- small_comps_by_net %>%
  summarise(
    n_redes_gt1     = dplyr::n(),
    n_small_total  = sum(n_small, na.rm = TRUE),
    mean_small_mu  = mean(mean_small, na.rm = TRUE),
    mean_small_sd  = sd(mean_small,   na.rm = TRUE),
    med_small_mu   = mean(med_small,  na.rm = TRUE),
    prop_isol_mu   = mean(prop_isol,  na.rm = TRUE),
    prop_isol_med  = median(prop_isol, na.rm = TRUE)
  )

print(small_comps_summary) # Resumen componentes pequeñas cuando n_components > 1

# hay componentes no aislados?
non_isolate_examples <- small_comps_by_net %>%
  filter(prop_isol < 1) %>%
  select(id, small_sizes, prop_isol) %>%
  head(10)

if (nrow(non_isolate_examples) > 0) {
  message("Ejemplos con pequeñas > 1 nodo:")
  print(non_isolate_examples)
} else {
  message("Todas las componentes pequeñas son aislados (tamaño 1).")
}

# Small-world por submuestreo (pendiente L ~ log N')
sw_scaling <- map_dfr(res_list, ~ .x$sw_sub %>% mutate(id = .x$id))
sw_fit <- sw_scaling %>% group_by(id) %>%
  do({
    df <- .
    mod <- lm(L ~ logN, data = df)
    tibble(slope = coef(mod)[["logN"]], r2 = summary(mod)$r.squared)
  }) %>% ungroup()

# Índice small-world sigma por red
sigma_tbl <- tibble(id = sapply(res_list, `[[`, "id"),
                    sigma = sapply(res_list, `[[`, "sigma"))

# Interpretación de sigma (Humphries & Gurney, 2008):
#   sigma = (C/C_ER) / (L/L_ER).
#   sigma > 1 indica estructura small-world: clustering mayor al de ER con longitudes de camino comparables.
#   En nuestras redes ATP, sigma ≈ 2.1 sugiere small-world robusto (alto C sin penalizar demasiado L).

# Rewiring bajo: deltas
rewiring_tbl <- map_dfr(res_list, ~{
  tibble(id = .x$id,
         L_rw = .x$rewiring["L_rw"],
         C_rw = .x$rewiring["C_rw"],
         rdeg_rw = .x$rewiring["rdeg_rw"])
}) %>% left_join(metrics_by_net %>% select(id, L= L, C_global, assortativity_degree = assort_deg),
                 by = "id") %>%
  mutate(
    dL = as.numeric(L_rw - L),
    dC = as.numeric(C_rw - C_global),
    dr = as.numeric(rdeg_rw - assortativity_degree)
  )

# Comparación de grados (ATP vs ER/WS/BA/CONFIG)
degree_cmp_tbl <- map_dfr(res_list, ~ .x$degree_cmp %>% mutate(id = .x$id))

# CCDF log-log “condensada” ATP vs BA
# (Se genera al vuelo; dejar plot comentado más abajo)

# Triadic census agregado (valores y z-score vs configuration model por red)
# Para z-score se necesitan nulos por red; arriba ya tenemos CONFIG como un modelo de referencia.
# Aquí calculamos proporciones por red y luego promediamos.
triad_by_net <- map_dfr(res_list, ~ .x$triad %>% mutate(id = .x$id))
triad_by_net <- triad_by_net %>%
  rowwise() %>%
  mutate(total = sum(c_across(starts_with("triads_"))),
         p0 = triads_0/total, p1 = triads_1/total, p2 = triads_2/total, p3 = triads_3/total) %>%
  ungroup()

triad_summary <- triad_by_net %>%
  summarise(across(c(p0,p1,p2,p3), list(mean=~mean(.x,na.rm=TRUE),
                                        median=~median(.x,na.rm=TRUE),
                                        sd=~sd(.x,na.rm=TRUE),
                                        p25=~quantile(.x,0.25,na.rm=TRUE),
                                        p75=~quantile(.x,0.75,na.rm=TRUE))))

# C(k) vs k condensado (mediana y bandas)
Ck_all <- map_dfr(res_list, ~ .x$Ck %>% mutate(id = .x$id))
Ck_summary <- Ck_all %>%
  group_by(k) %>%
  summarise(Ck_med = median(Ck, na.rm = TRUE),
            Ck_p25 = quantile(Ck, 0.25, na.rm = TRUE),
            Ck_p75 = quantile(Ck, 0.75, na.rm = TRUE),
            n = sum(!is.na(Ck)),
            .groups = "drop")

# Percolación condensada
percolation_all <- bind_rows(lapply(res_list, `[[`, "percolation"))
percolation_summary <- percolation_all %>%
  group_by(mode, p) %>%
  summarise(S_med = median(S, na.rm = TRUE),
            S_p25 = quantile(S, 0.25, na.rm = TRUE),
            S_p75 = quantile(S, 0.75, na.rm = TRUE),
            .groups = "drop")

# ---------- Escritura de tablas (comentadas) ----------
# dir.create("results", showWarnings = FALSE)
# write_csv(metrics_by_net,   "results/metrics_by_network.csv")
# write_csv(metrics_summary,  "results/metrics_summary.csv")
# write_csv(sw_scaling,       "results/smallworld_subsampling.csv")
# write_csv(sw_fit,           "results/smallworld_subsampling_fits.csv")
# write_csv(sigma_tbl,        "results/smallworld_sigma.csv")
# write_csv(rewiring_tbl,     "results/rewiring_effects.csv")
# write_csv(degree_cmp_tbl,   "results/degree_comparisons.csv")
# write_csv(triad_by_net,     "results/triadic_census_by_network.csv")
# write_csv(triad_summary,    "results/triadic_census_summary.csv")
# write_csv(Ck_all,           "results/Ck_by_network.csv")
# write_csv(Ck_summary,       "results/Ck_summary.csv")
# write_csv(percolation_all,  "results/percolation_by_network.csv")
# write_csv(percolation_summary, "results/percolation_summary.csv")


# ==============================================================================
# PLOTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 0) Distribuciones de grado
# ------------------------------------------------------------------------------

# Colores solicitados
col_ER <- "#2ECC71"   # verde
col_BA <- "#6A1B9A"   # morado
col_ATP <- "#E74C3C"  # rojo

# Para ATP: condensado sobre 100 redes (media y sd de p_k por k).
# Para ER (G(N,p)): Binomial(N-1, p) con p = <k>/(N-1) = 2m / [N(N-1)], usando m0 = round(mean(m))
# Para BA: P(k) = 2 m (m+1) / [k (k+1) (k+2)] para k >= m; 0 en otro caso.
N0 <- metrics_by_net$N[1]
m0 <- round(mean(metrics_by_net$m[1:100]))
k_mean0 <- 2 * m0 / N0
p_er <- (2 * m0) / (N0 * (N0 - 1))  # => <k> = p*(N-1) ≈ 2m/N

## Condensar ATP: media y sd por k (promediando sobre las 100 redes)
## Nota: completamos los grados faltantes por red con cero para no sesgar la media.
deg_atp_df <- map_dfr(res_list, ~{
  tibble(k = as.integer(degree(.x$g)))
}, .id = "id")

k_max_atp <- max(deg_atp_df$k, na.rm = TRUE)

# Conteo por red y grado
deg_atp_counts <- deg_atp_df %>%
  group_by(id, k) %>%
  summarise(n = n(), .groups = "drop")

# Completar grados faltantes 0..k_max_atp con n=0 por red
deg_atp_counts_full <- deg_atp_counts %>%
  group_by(id) %>%
  tidyr::complete(k = 0:k_max_atp, fill = list(n = 0L)) %>%
  ungroup()

# pk por red
deg_atp_pk <- deg_atp_counts_full %>%
  group_by(id) %>%
  mutate(pk = n / sum(n)) %>%
  ungroup()

# Media y desvío estándar por k (sobre redes)
deg_atp_pk_stats <- deg_atp_pk %>%
  group_by(k) %>%
  summarise(mean_pk = mean(pk, na.rm = TRUE),
            sd_pk   = sd(pk,   na.rm = TRUE),
            .groups = "drop")

# ER teórico: Binomial(N0-1, p_er)
k_er <- 0:(N0 - 1)
er_theory <- tibble(
  model = "ER (teórico)",
  k = k_er,
  pk = dbinom(k_er, size = N0 - 1, prob = p_er)
)

# BA teórico: P(k) = 2 m (m+1) / [k (k+1) (k+2)] para k >= m, 0 otro caso

m_ba <- max(1L, round(k_mean0 / 2))           # m_BA ~ <k>/2
k_ba <- m_ba:(max(k_er))                      # soporte desde m hasta N-1
ba_theory <- tibble(
  model = "BA (teórico)",
  k = k_ba,
  pk = (2 * m_ba * (m_ba + 1)) / (k_ba * (k_ba + 1) * (k_ba + 2))
)

# Guía power-law (solo para k < 14), gamma = 2.81, normalizada en k0 = 14
gamma_left <- 2.81
k0_left <- 14L
# Valor BA exacto en k0 (si m_ba > k0, usamos m_ba para evitar valores 0; pero por pedido, forzamos k0=14)
P_ba_k0 <- if (k0_left >= m_ba) {
  (2 * m_ba * (m_ba + 1)) / (k0_left * (k0_left + 1) * (k0_left + 2))
} else {
  # Si k0 < m_ba, BA exacta vale 0; en ese caso, normalizamos usando k = m_ba
  (2 * m_ba * (m_ba + 1)) / (m_ba * (m_ba + 1) * (m_ba + 2))
}
C_left <- P_ba_k0 * (k0_left ^ gamma_left)
k_left <- 1:(k0_left - 1L)
ba_power_left <- tibble(
  k  = k_left,
  pk = C_left * (k_left ^ (-gamma_left))
) %>% filter(pk > 0)

# (a) Escala lineal
p_lin <- ggplot() +
  geom_line(data = ba_theory, aes(x = k, y = pk, color = "BA (teórico)"), linewidth = 0.9, show.legend = FALSE) +
  geom_line(data = er_theory, aes(x = k, y = pk, color = "ER (teórico)"), linewidth = 0.9, show.legend = FALSE) +
  geom_line(data = deg_atp_pk_stats, aes(x = k, y = mean_pk, color = "ATP (media)"), linewidth = 1, show.legend = FALSE) +
  geom_ribbon(data = deg_atp_pk_stats,
              aes(x = k, ymin = pmax(mean_pk - sd_pk, 0), ymax = mean_pk + sd_pk, fill = "ATP (±1 sd)"),
              alpha = 0.5, inherit.aes = FALSE, show.legend = FALSE) +
  geom_line(data = ba_power_left, aes(x = k, y = pk),
            color = col_BA, linetype = "dashed", linewidth = 0.9) +
  scale_color_manual(values = c("BA (teórico)" = col_BA,
                                "ER (teórico)" = col_ER,
                                "ATP (media)"  = col_ATP)) +
  scale_fill_manual(values  = c("ATP (±1 sd)" = col_ATP)) +
  labs(title = "(a) Escala lineal",
       x = "k", y = "p_k", color = NULL, fill = NULL) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.08)) +
  theme_minimal()

# (b) Escala log–log
er_theory_log <- er_theory %>% filter(k >= 1, pk > 0)
ba_theory_log <- ba_theory %>% filter(k >= 1, pk > 0)
atp_log <- deg_atp_pk_stats %>% filter(k >= 1, mean_pk > 0)

p_log <- ggplot() +
  geom_line(data = ba_theory_log, aes(x = k, y = pk, color = "BA (teórico)"), linewidth = 0.9) +
  geom_line(data = er_theory_log, aes(x = k, y = pk, color = "ER (teórico)"), linewidth = 0.9) +
  geom_line(data = atp_log, aes(x = k, y = mean_pk, color = "ATP (media)"), linewidth = 1) +
  geom_ribbon(data = atp_log,
              aes(x = k, ymin = pmax(mean_pk - sd_pk, 0), ymax = mean_pk + sd_pk, fill = "ATP (±1 sd)"),
              alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = ba_power_left %>% filter(k >= 1, pk > 0),
            aes(x = k, y = pk),
            color = col_BA, linetype = "dashed", linewidth = 0.9) +
  scale_x_log10() + scale_y_log10() +
  coord_cartesian(xlim = c(5, 150), ylim = c(1e-5, 1)) +
  scale_color_manual(values = c("BA (teórico)" = col_BA,
                                "ER (teórico)" = col_ER,
                                "ATP (media)"  = col_ATP)) +
  scale_fill_manual(values  = c("ATP (±1 sd)" = col_ATP)) +
  labs(title = "(b) Escala log–log",
       x = "k (log)", y = "p_k (log)", color = NULL, fill = NULL) +
  theme_minimal()

# Combinar en una sola figura (como la referencia adjunta)
(p_lin | p_log) #+ plot_annotation(title = "Distribución de grado: ATP (rojo) vs ER (verde) y BA (morado; teóricos)")
ggsave("trabajo_1_plots/ATP_network_tests/degree_dist.pdf", width = 11, height = 4.8, dpi = 300)

# ------------------------------------------------------------------------------
# 1) Propiedad Small-world
# ------------------------------------------------------------------------------

# A) Small-world por submuestreo (condensado)
# Deducción: L_ER ≈ log(N)/log(<k>). Con N=1000 y <k>≈28 → L_ER ≈ log(1000)/log(28) ≈ 2.6
L_ER <- 2.6
ggplot(sw_scaling, aes(x = logN, y = L, group = id)) +
  geom_line(alpha = 0.15) +
  geom_smooth(se = TRUE, color = "black") +
  geom_hline(yintercept = L_ER, linetype = "dotted", color = "red") +
  annotate("text", x = min(sw_scaling$logN), y = L_ER, vjust = -0.5, hjust = 0,
           label = "Baseline ER: L_ER ≈ 2.6") +
  labs(title = "Small-world por submuestreo (L vs log N')",
       x = "log(N')", y = "Longitud media de camino (L)") +
  theme_minimal()
# ggsave("trabajo_1_plots/ATP_network_tests/smallworld_subsampling.pdf", width = 7, height = 4.5, dpi = 300)

#
# B) Rewiring: antes vs después
p1 <- metrics_by_net %>%
  select(id, L, C_global, assort_deg) %>%
  left_join(rewiring_tbl %>% select(id, L_rw, C_rw, rdeg_rw), by = "id")
ggplot(filter(p1, is.finite(L) & is.finite(L_rw))) +
  geom_point(aes(x = L, y = L_rw), alpha = 0.5) +
  geom_abline(linetype = "dashed") +
  labs(title = "Efecto del rewiring bajo en L", x = "L (antes)", y = "L (después)") +
  theme_minimal()
# ggsave("trabajo_1_plots/ATP_network_tests/rewiring_L.pdf", width = 5, height = 4, dpi = 300)

ggplot(p1) +
  geom_point(aes(x = C_global, y = C_rw), alpha = 0.5) +
  geom_abline(linetype = "dashed") +
  labs(title = "Efecto del rewiring bajo en C_global", x = "C antes", y = "C después") +
  theme_minimal()
# ggsave("trabajo_1_plots/ATP_network_tests/rewiring_C.pdf", width = 5, height = 4, dpi = 300)

# ------------------------------------------------------------------------------
# 2) ATP vs BA
# ------------------------------------------------------------------------------

# C) CCDF log-log de grados: ATP (condensado) vs BA
deg_df <- map_dfr(res_list, ~{
  d <- degree(.x$g)
  tibble(k = as.integer(d))
}, .id = "id")

# CCDF condensada (media por k). Para evitar ruido en k muy bajos, filtramos k>=1
ccdf_atp <- deg_df %>%
  filter(k >= 1) %>%
  group_by(k) %>%
  summarise(CCDF = mean(1 - ecdf(k)(k), na.rm = TRUE), .groups="drop") %>%
  mutate(model = "ATP (100 redes)")

# BA de referencia usando N y m de la primera red
N0 <- metrics_by_net$N[1]; m0 <- metrics_by_net$m[1]
g_ba0 <- ba_match(N0, m0)
deg_ba <- degree(g_ba0)
ba_tbl <- tibble(k = sort(unique(as.integer(deg_ba)))) %>%
  mutate(CCDF = 1 - ecdf(deg_ba)(k),
         model = "Barabási–Albert")

ccdf_plot_df <- bind_rows(ccdf_atp, ba_tbl)

ggplot(ccdf_plot_df, aes(x = k, y = CCDF, color = model)) +
  geom_point(data = subset(ccdf_plot_df, model == "ATP (100 redes)"), alpha = 0.6) +
  geom_line(data = subset(ccdf_plot_df, model == "Barabási–Albert"), linewidth = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "CCDF log-log: ATP (condensada) vs BA",
       x = "k (log)", y = "P(K≥k) (log)", color = NULL) +
  theme_minimal()
# ggsave("trabajo_1_plots/ATP_network_tests/degree_ccdf_loglog_ATP_vs_BA.pdf", width = 6.5, height = 4.5, dpi = 300)

# ---------- Test formal de cola larga (power-law) ----------
#   - p_value grande (>= 0.1) -> no podemos rechazar power-law (plausible).
#   - p_value pequeño (< 0.1) -> rechazamos power-law (no plausible); considerar lognormal/Weibull.

library(poweRlaw)
deg_all <- unlist(lapply(res_list, function(x) degree(x$g)), use.names = FALSE)
deg_all <- deg_all[deg_all > 0]  # poweRlaw (discreto) requiere valores positivos
m_pl <- displ$new(deg_all)
est  <- estimate_xmin(m_pl)
m_pl$setXmin(est)

# KS bootstrap para p-value
bs   <- bootstrap_p(m_pl, no_of_sims = 500, threads = 8)
heavy_tail_results <- list(
  xmin  = m_pl$getXmin(),
  alpha = m_pl$pars,
  p_value = bs$p
)
print(heavy_tail_results)
# Comparación vs lognormal (opcional)
m_ln <- dislnorm$new(deg_all)
m_ln$setXmin(m_pl$getXmin())
cmp <- compare_distributions(m_pl, m_ln)
print(cmp)

# ------------------------------------------------------------------------------
# 3) Censo triádico
# ------------------------------------------------------------------------------

# D) Censo triádico (sin p0) + referencia ER
# En ER con probabilidad de arista p ≈ densidad, las probabilidades por triada son:
#   P(3 aristas) = p^3
#   P(2 aristas) = 3 p^2 (1-p)
#   P(1 arista)  = 3 p (1-p)^2
# Usamos p_bar = mean(density) como baseline global.
p_bar <- mean(metrics_by_net$density, na.rm = TRUE)
triad_refs <- tibble(
  triad = c("p1","p2","p3"),
  ref   = c(3 * p_bar * (1 - p_bar)^2,
            3 * p_bar^2 * (1 - p_bar),
            p_bar^3)
)
triad_long <- triad_by_net %>%
  transmute(id, p1, p2, p3) %>%
  pivot_longer(-id, names_to = "triad", values_to = "prop")

# ---- Resumen a consola para comparar con ER ----
triad_medians <- triad_long %>%
  group_by(triad) %>%
  summarise(median_prop = median(prop, na.rm = TRUE),
            p25 = quantile(prop, 0.25, na.rm = TRUE),
            p75 = quantile(prop, 0.75, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(triad_refs, by = "triad") %>%
  mutate(delta_vs_ER = median_prop - ref,
         ratio_vs_ER = median_prop / ref)
print(triad_medians)
# En redes sociales reales, se espera:
# - Exceso de triángulos (triads_3) y de wedges (triads_2) respecto a grafos aleatorios con igual N y m,
#   reflejando cierre triádico / clustering no trivial.
# - Déficit de triadas con 0 aristas (triads_0) frente a nulos.
# - Este patrón debería mantenerse incluso controlando por la secuencia de grados (configuration model),
#   donde el exceso de triángulos suele seguir siendo positivo si hay homofilia y geometría social.
# (Síntesis basada en Talaga & Nowak, 2020).
# -----------------------------------------------

ggplot(triad_long, aes(x = triad, y = prop)) +
  stat_summary(fun = median, geom = "col", fill = "grey70") +
  stat_summary(fun.min = function(z) quantile(z,0.25), fun.max = function(z) quantile(z,0.75),
               geom = "errorbar", width = 0.3) +
  geom_hline(data = triad_refs, aes(yintercept = ref), color = "red", linetype = "dashed") +
  labs(title = "Censo triádico (ATP): p1, p2, p3 con referencia ER (línea roja)",
       x = "Tipo de triada", y = "Proporción") +
  theme_minimal()
# Referencia ER: p_bar = mean(density); p1=3p(1-p)^2, p2=3p^2(1-p), p3=p^3.
# ggsave("trabajo_1_plots/ATP_network_tests/triadic_census_summary.pdf", width = 7, height = 4, dpi = 300)

# ------------------------------------------------------------------------------
# 4) C(k) vs k y percolación
# ------------------------------------------------------------------------------

# E) C(k) vs k (condensado)
# Deducción: en ER, C ≈ p y p ≈ <k>/N → C_ER ≈ 28/1000 = 0.028
C_ER <- 28/1000
ggplot(Ck_summary, aes(x = k, y = Ck_med)) +
  geom_ribbon(aes(ymin = Ck_p25, ymax = Ck_p75), fill = "#E74C3C", alpha = 0.5, color = NA) +
  geom_line(color = "#E74C3C", linewidth = 1.1) +
  geom_hline(yintercept = C_ER, linetype = "dashed", color = "#2ECC71", linewidth = 1.1) +
  annotate("text", x = min(Ck_summary$k, na.rm = TRUE), y = C_ER, vjust = -0.51, hjust = -1.5,
           label = "C_ER = <k>/N = 0.028", color = "#2ECC71") +
  labs(title = "C(k) vs k (ATP, mediana y bandas [p25,p75])",
       x = "k", y = "C(k)") + theme_minimal()
ggsave("trabajo_1_plots/ATP_network_tests/Ck_vs_k.pdf", width = 6.5, height = 4, dpi = 300)

# print(Ck_summary %>% arrange(desc(k)) %>% select(k, n) %>% head(10)) # conteo colas

# F) Percolación: fallas aleatorias vs ataques dirigidos (condensado)
ggplot(percolation_summary, aes(x = p, y = S_med, color = mode)) +
  geom_line(linewidth = 1) +
  labs(title = "Percolación del GCC: fallas aleatorias vs ataques por grado",
       x = "Fracción removida (p)", y = "Tamaño relativo GCC") +
  theme_minimal()
# ggsave("trabajo_1_plots/ATP_network_tests/percolation_curves.pdf", width = 6.5, height = 4, dpi = 300)



# ==============================================================================
# REPORTE
# ==============================================================================

# ---------- Utilidades de formateo ----------
fmt_pm <- function(mu, sd, digits = 3) {
  if (is.na(mu) || is.na(sd)) return("--")
  paste0(format(round(mu, digits), nsmall = digits), " $\\pm$ ", format(round(sd, digits), nsmall = digits))
}
fmt_num <- function(x, digits = 3) {
  if (is.na(x)) return("--")
  format(round(x, digits), nsmall = digits)
}

# ---------- Rescatar/derivar estadísticas ----------
# metrics_summary ya contiene mean/median/sd/p25/p75 por métrica
ms <- metrics_summary %>% select(metric, mean, sd, median, p25, p75) %>% tidyr::pivot_wider(names_from = metric, values_from = c(mean, sd, median, p25, p75))

# Sigma (small-world): media y sd
sigma_mu <- mean(sigma_tbl$sigma, na.rm = TRUE)
sigma_sd <- sd(sigma_tbl$sigma, na.rm = TRUE)

# Triadic census: medianas por p1,p2,p3 (sobre 100 redes)
tri_mu_med <- triad_summary %>% transmute(
  p1_med = p1_median, p2_med = p2_median, p3_med = p3_median,
  p1_sd  = p1_sd,     p2_sd  = p2_sd,     p3_sd  = p3_sd
)
p1_med <- tri_mu_med$p1_med; p2_med <- tri_mu_med$p2_med; p3_med <- tri_mu_med$p3_med
p1_sd  <- tri_mu_med$p1_sd;  p2_sd  <- tri_mu_med$p2_sd;  p3_sd  <- tri_mu_med$p3_sd

# Referencias ER (basadas en densidad promedio p_bar)
p_bar <- mean(metrics_by_net$density, na.rm = TRUE)
p1_ER <- 3 * p_bar * (1 - p_bar)^2
p2_ER <- 3 * p_bar^2 * (1 - p_bar)
p3_ER <- p_bar^3

# Ratios triádicos vs ER (usar medianas ATP)
r1 <- p1_med / p1_ER
r2 <- p2_med / p2_ER
r3 <- p3_med / p3_ER

# Baselines ER adicionales
N0_rep  <- as.integer(ms$median_N)
kbar_rep <- ms$mean_k_mean
dens_rep <- ms$mean_density
L_ER_baseline <- log(N0_rep) / log(kbar_rep)  # aproximación clásica L_ER ≈ log N / log <k>
C_ER_baseline <- kbar_rep / N0_rep           # C_ER ≈ <k> / N
# Nota: estos baselines se informan al pie de la tabla.

gini_summary <- degree_cmp_tbl %>% 
  filter(model %in% c("ER","BA")) %>% 
  group_by(model) %>% 
  summarise(
    mean_gini   = mean(gini, na.rm = TRUE),
    sd_gini     = sd(gini, na.rm = TRUE),
    median_gini = median(gini, na.rm = TRUE),
    .groups = "drop"
  )

gini_BA <- gini_summary$mean_gini[1]
gini_ER <- gini_summary$mean_gini[2]

# ---------- Construir LaTeX (standalone) ----------
tex_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=2.5cm]{geometry}",
  "\\usepackage{booktabs}",
  "\\usepackage{siunitx}",
  "\\usepackage{amsmath}",
  "\\usepackage{array}",
  "\\usepackage{caption}",
  "\\captionsetup[table]{labelfont=bf}",
  "\\begin{document}",
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Resumen topológico de 100 redes ATP (N=1000). Medias y desviaciones estándar entre réplicas. Se excluyen comparaciones de distribución de grados por requerimiento.}",
  "\\vspace{0.5em}",
  "\\begin{tabular}{l@{\\hspace{1em}}c@{\\hspace{3em}}l}",
  "\\toprule",
  "\\textbf{Métrica} & \\textbf{Valor resumen (ATP)} & \\textbf{Referencia / Interpretación} \\\\",
  "\\midrule",
  sprintf("N & %s & tamaño de red (fijo). \\\\", fmt_num(ms$median_N, 0)),
  sprintf("$\\langle k \\rangle$ & %s & grado medio. \\\\", fmt_pm(ms$mean_k_mean, ms$sd_k_mean, 2)),
  sprintf("Densidad & %s & fracción de enlaces presentes. \\\\", fmt_pm(ms$mean_density, ms$sd_density, 4)),
  sprintf("$L$ (camino medio) & %s & pequeño; cf. $L_{\\mathrm{ER}}\\approx$ %s. \\\\", fmt_pm(ms$mean_L, ms$sd_L, 3), fmt_num(L_ER_baseline, 2)),
  sprintf("$D$ (diámetro) & %s & muy bajo. \\\\", fmt_pm(ms$mean_D, ms$sd_D, 2)),
  sprintf("$C_{\\mathrm{global}}$ & %s & $\\gg C_{\\mathrm{ER}}\\approx$ %s. \\\\", fmt_pm(ms$mean_C_global, ms$sd_C_global, 4), fmt_num(C_ER_baseline, 3)),
  sprintf("$C_{\\mathrm{avg}}$ & %s & coherente con clustering elevado. \\\\", fmt_pm(ms$mean_C_avg, ms$sd_C_avg, 4)),
  sprintf("Asortatividad (grado) & %s & positiva (homofilia por grado). \\\\", fmt_pm(ms$mean_assort_deg, ms$sd_assort_deg, 3)),
  sprintf("$G$ (Gini grado) & %s & Heterog. moderada; cf. $G_{\\mathrm{ER}}\\approx %s$, $G_{\\mathrm{SF}}\\approx %s$. \\\\", 
          fmt_pm(ms$mean_gini_deg, ms$sd_gini_deg, 3),
          fmt_num(gini_ER, 2),
          fmt_num(gini_BA, 2)),
  sprintf("Small-world $\\sigma$ & %s & $(C/C_{\\mathrm{ER}})/(L/L_{\\mathrm{ER}})$. \\\\", fmt_pm(sigma_mu, sigma_sd, 2)),
  "\\midrule",
  sprintf("$p_1$ (1 arista) & %s & ER: %s; $p_1/p_1^{\\mathrm{ER}}$ = %s. \\\\",
          fmt_pm(p1_med, p1_sd, 4), fmt_num(p1_ER, 4), fmt_num(r1, 2)),
  sprintf("$p_2$ (2 aristas) & %s & ER: %s; $p_2/p_2^{\\mathrm{ER}}$ = %s. \\\\",
          fmt_pm(p2_med, p2_sd, 5), fmt_num(p2_ER, 5), fmt_num(r2, 2)),
  sprintf("$p_3$ (3 aristas) & %s & ER: %s; $p_3/p_3^{\\mathrm{ER}}$ = %s. \\\\",
          fmt_pm(p3_med, p3_sd, 6), fmt_num(p3_ER, 6), fmt_num(r3, 2)),
  "\\end{tabular}",
  "\\vspace{0.75em}",
  "\\captionsetup{justification=raggedright,singlelinecheck=false}",
  "\\small \\textit{Notas de referencia y cálculos:}",
  "\\begin{itemize}",
  sprintf("\\item $p=\\overline{\\text{density}}= %s$. Para ER: $p_1=3p(1-p)^2$, $p_2=3p^2(1-p)$, $p_3=p^3$.", fmt_num(p_bar, 4)),
  sprintf("\\item $C_{\\mathrm{ER}}\\approx \\langle k \\rangle/N = %s / %s = %s$ (aprox.).", fmt_num(kbar_rep, 2), fmt_num(N0_rep, 0), fmt_num(C_ER_baseline, 3)),
  sprintf("\\item $L_{\\mathrm{ER}}\\approx \\log N / \\log \\langle k \\rangle = \\log(%s)/\\log(%s) \\approx %s$ (aprox.).", fmt_num(N0_rep, 0), fmt_num(kbar_rep, 2), fmt_num(L_ER_baseline, 2)),
  "\\item $\\sigma = (C/C_{\\mathrm{ER}})/(L/L_{\\mathrm{ER}})$ (Humphries \\& Gurney, 2008).",
  "\\end{itemize}",
  "\\end{table}",
  "\\end{document}"
)

# Escribir archivo .tex
writeLines(tex_lines, con = "trabajo_1_files/ATP_network_tests/tabla_resumen_ATP.tex")

message('Archivo LaTeX escrito en: results/tabla_resumen_ATP.tex. Compila con: pdflatex results/tabla_resumen_ATP.tex')
