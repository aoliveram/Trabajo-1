## ==== Regressions (final): GRAND_COMBINED-only, robust load, BAM ==========

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr); library(stringr); library(mgcv)
})

# --------------------------- 0) Parámetros ---------------------------------
NUM_CORES_TO_USE <- 8L
K_BIVAR          <- 10L        # base dimension por superficie 2D
OUT_DIR_MODELS   <- "models"
OUT_DIR_FIGS     <- "figs"
dir.create(OUT_DIR_MODELS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_FIGS,   showWarnings = FALSE, recursive = TRUE)

# --------------------------- 1) Directorios --------------------------------
dir_atp_emp <- "trabajo_1_files/ATP_diffusion_simulation_files_sigm"
dir_atp_er  <- "trabajo_1_files/ATP_ER_diffusion_simulation_files_sigm"
dir_gss_emp <- "trabajo_1_files/GSS_diffusion_simulation_files_sigm"
dir_gss_er  <- "trabajo_1_files/GSS_ER_diffusion_simulation_files_sigm"

ls_grand_combined <- function(path) {
  if (!dir.exists(path)) return(character(0))
  f <- list.files(path, pattern = "\\.rds$", full.names = TRUE, recursive = FALSE)
  f[grepl("GRAND_COMBINED", basename(f), ignore.case = TRUE)]
}

files_atp_emp <- ls_grand_combined(dir_atp_emp)
files_atp_er  <- ls_grand_combined(dir_atp_er)
files_gss_emp <- ls_grand_combined(dir_gss_emp)
files_gss_er  <- ls_grand_combined(dir_gss_er)

stopifnot(length(files_atp_emp) + length(files_atp_er) +
            length(files_gss_emp) + length(files_gss_er) > 0)

make_sources <- function(paths, contagion_type, network_type) {
  tibble(path = paths,
         contagion_type = contagion_type,
         network_type   = network_type)
}

sources <- bind_rows(
  make_sources(files_atp_emp, "Innovation",       "Empirical"),
  make_sources(files_atp_er,  "Innovation",       "ER"),
  make_sources(files_gss_emp, "CollectiveAction", "Empirical"),
  make_sources(files_gss_er,  "CollectiveAction", "ER")
)

message("== Discovered GRAND_COMBINED files ==")
print(sources, n = 50)

# -------------------- 2) Helpers: parse / extract --------------------------
parse_key_into_cols <- function(df, key){
  if (is.null(key) || is.na(key) || !nzchar(key)) return(df)
  m_sd <- str_match(key, "sd[_]?([0-9.]+)")
  if (!is.na(m_sd[1,2])) df$tau_sd <- dplyr::coalesce(df$tau_sd, as.numeric(m_sd[1,2]))
  m_mu <- str_match(key, "(mean|mu)[_]?([0-9.]+)")
  if (!is.na(m_mu[1,3])) df$tau_mean <- dplyr::coalesce(df$tau_mean, as.numeric(m_mu[1,3]))
  df
}

collect_dfs <- function(x, key_path = character()){
  out <- list()
  if (inherits(x, "data.frame")) {
    df <- x
    if (length(key_path)) {
      df$.nest_path <- paste(key_path, collapse = "/")
      for (k in key_path) df <- parse_key_into_cols(df, k)
    }
    out <- list(df)
  } else if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) {
      for (i in seq_along(x)) out <- c(out, collect_dfs(x[[i]], c(key_path, paste0("idx", i))))
    } else {
      for (nm in nms) out <- c(out, collect_dfs(x[[nm]], c(key_path, nm)))
    }
  }
  out
}

parse_seed_from_filename <- function(path) {
  b <- basename(path)
  s1 <- str_match(b, "means_([A-Za-z]+)\\.rds$")[,2]
  s2 <- str_match(b, "_(random|central|eigen|closeness|marginal)\\.rds$")[,2]
  seed <- dplyr::coalesce(s1, s2)
  if (is.na(seed)) seed <- "random"
  seed
}

# ---------------- 3) Lectura robusta + armonización ------------------------
safe_read_and_label <- function(path, contagion_type, network_type) {
  stopifnot(file.exists(path))
  obj <- readRDS(path)
  dfs <- collect_dfs(obj)
  if (!length(dfs)) stop("No data.frames in RDS: ", path)
  
  out <- suppressMessages(bind_rows(dfs))
  
  # Renombres canónicos
  if ("innovation_iul_Gamma" %in% names(out)) out <- rename(out, IUL = innovation_iul_Gamma)
  if ("social_distance_h"    %in% names(out)) out <- rename(out, MSP = social_distance_h)
  
  # Proporciones (si solo hay conteos)
  if (!"adopt_total" %in% names(out) && all(c("num_adopters","N_nodes_actual") %in% names(out))) {
    out$adopt_total <- out$num_adopters / out$N_nodes_actual
  }
  if (!"adopt_social" %in% names(out) && all(c("num_adopted_social","N_nodes_actual") %in% names(out))) {
    out$adopt_social <- out$num_adopted_social / out$N_nodes_actual
  }
  if (!"adopt_rational" %in% names(out) && all(c("num_adopted_rational","N_nodes_actual") %in% names(out))) {
    out$adopt_rational <- out$num_adopted_rational / out$N_nodes_actual
  }
  
  # Thresholds si llegan con otros nombres
  if (!"tau_mean" %in% names(out) && "threshold_mean_param" %in% names(out)) out$tau_mean <- out$threshold_mean_param
  if (!"tau_sd"   %in% names(out) && "threshold_sd_param"   %in% names(out)) out$tau_sd   <- out$threshold_sd_param
  
  # seed_type si falta
  if (!"seed_type" %in% names(out)) out$seed_type <- parse_seed_from_filename(path)
  
  out$contagion_type <- contagion_type
  out$network_type   <- network_type
  out
}

df_all <- purrr::pmap_dfr(sources, ~ safe_read_and_label(..1, ..2, ..3))

# Tipado y limpieza
df_all <- df_all %>%
  mutate(
    IUL            = as.numeric(IUL),
    MSP            = as.numeric(MSP),
    tau_mean       = as.numeric(tau_mean),
    tau_sd         = as.numeric(tau_sd),
    adopt_total    = as.numeric(adopt_total),
    adopt_social   = suppressWarnings(as.numeric(adopt_social)),
    adopt_rational = suppressWarnings(as.numeric(adopt_rational)),
    seed_type      = as.factor(seed_type),
    contagion_type = factor(contagion_type, levels = c("Innovation","CollectiveAction")),
    network_type   = factor(network_type,   levels = c("Empirical","ER"))
  ) %>%
  filter(is.finite(IUL), is.finite(MSP),
         is.finite(tau_mean), is.finite(tau_sd),
         is.finite(adopt_total),
         adopt_total >= 0, adopt_total <= 1)

# Baseline: random
df_all <- df_all %>% mutate(seed_type = relevel(seed_type, ref = "random"))

# Cobertura
message("\n=== Coverage check (GRAND_COMBINED only) ===")
message("N rows: ", nrow(df_all))
message("\ncontagion_type:\n"); print(table(df_all$contagion_type, useNA = "ifany"))
message("\nnetwork_type:\n");   print(table(df_all$network_type,   useNA = "ifany"))
message("\nseed_type:\n");      print(table(df_all$seed_type,      useNA = "ifany"))

has_multi_seed       <- nlevels(df_all$seed_type)      >= 2
has_both_contagions  <- nlevels(df_all$contagion_type) >= 2
has_both_networks    <- nlevels(df_all$network_type)   >= 2

# ---------------- 4) Construcción segura de fórmulas -----------------------
smooth_term <- if (has_both_contagions) {
  paste0("s(IUL, MSP, k=", K_BIVAR, ", by=contagion_type)")
} else {
  paste0("s(IUL, MSP, k=", K_BIVAR, ")")
}

rhs_common <- c(
  smooth_term,
  "tau_mean + tau_sd",
  if (has_multi_seed) "seed_type" else NULL,
  if (has_both_networks) "network_type" else NULL
)

form_total  <- as.formula(paste("adopt_total  ~", paste(rhs_common, collapse = " + ")))
form_social <- as.formula(paste("adopt_social ~", paste(rhs_common, collapse = " + ")))
form_rational <- as.formula(paste("adopt_rational ~", paste(rhs_common, collapse = " + ")))

message("Formula TOTAL : ", deparse(form_total))
message("Formula SOCIAL: ", deparse(form_social))
message("Formula RATIONAL: ", deparse(form_rational))

# ---------------- 5) Función de ajuste BAM + guardado ----------------------
fit_bam <- function(formula, data, label, out_prefix){
  ctrl <- gam.control(nthreads = NUM_CORES_TO_USE)
  t0 <- Sys.time()
  m  <- bam(
    formula  = formula,
    family   = binomial(link="logit"),
    data     = data,
    method   = "fREML",
    discrete = TRUE,
    nthreads = NUM_CORES_TO_USE
  )
  t1 <- Sys.time()
  elapsed <- difftime(t1, t0, units = "secs")
  msg <- paste0("[", label, "] elapsed (s): ", round(as.numeric(elapsed), 2))
  message(msg)
  
  # Guardar modelo y resumen
  saveRDS(m, file = file.path(OUT_DIR_MODELS, paste0(out_prefix, ".rds")))
  sink(file.path(OUT_DIR_MODELS, paste0(out_prefix, "_summary.txt")))
  on.exit(sink(), add = TRUE)
  cat("\n== ", label, " ==\n")
  print(summary(m)$p.table)     # solo los coeficientes paramétricos
  print(summary(m)$s.table)     # solo resumen de smooths
  cat("\n-- Deviance explained:", summary(m)$dev.expl, "\n")
  sink()  # cerrar
}

fit_bam <- function(formula, data, label, out_prefix){
  ctrl <- gam.control(nthreads = NUM_CORES_TO_USE)
  t0 <- Sys.time()
  m  <- bam(
    formula  = formula,
    family   = binomial(link="logit"),
    data     = data,
    method   = "fREML",
    discrete = TRUE,
    nthreads = NUM_CORES_TO_USE
    # Si mantienes proporciones y quieres quitar el warning, añade:
    # , weights = if ("N_nodes_actual" %in% names(data)) data$N_nodes_actual else NULL
  )
  t1 <- Sys.time()
  elapsed <- round(as.numeric(difftime(t1, t0, units = "secs")), 2)
  message("[", label, "] elapsed (s): ", elapsed)
  
  # 1) Mostrar summary COMPLETO en consola
  sm <- summary(m)
  cat("\n== ", label, " ==\n", sep = "")
  print(sm)
  
  # 2) Guardar un resumen compacto a archivo (paramétricos + smooths + dev.expl)
  out_txt <- file.path(OUT_DIR_MODELS, paste0(out_prefix, "_summary.txt"))
  capture.output({
    cat("\n== ", label, " ==\n", sep = "")
    print(sm$p.table)  # coeficientes paramétricos
    cat("\n-- smooths --\n")
    print(sm$s.table)  # resumen de suavizados
    cat("\nDeviance explained: ", sm$dev.expl, "\n", sep = "")
    cat("Elapsed (s): ", elapsed, "\n", sep = "")
  }, file = out_txt, type = "output")
  
  # 3) Guardar el modelo
  saveRDS(m, file.path(OUT_DIR_MODELS, paste0(out_prefix, ".rds")))
  
  invisible(m)
}

# ---------------- 6) Ajustes FAMILIA A (Total) y B (Social) ----------------
m_bam_total <- fit_bam(form_total,  df_all, "BAM Total Adoption",  "bam_total")
m_bam_social <- fit_bam(form_social, df_all, "BAM Social Adoption", "bam_social")
m_bam_rational <- fit_bam(form_rational, df_all, "BAM Rational Adoption", "bam_rational")

# ---------------- 7) Figuras (superficies 2D) ------------------------------
# Guardamos PNGs para no depender de dispositivo gráfico interactivo
plot_surface <- function(model, fname_prefix){
  png(file.path(OUT_DIR_FIGS, paste0(fname_prefix, "_surfaces.png")),
      width = 1600, height = 800, res = 150)
  plot(model, scheme = 2, pages = 1, shade = TRUE)
  dev.off()
}
plot_surface(m_bam_total,  "bam_total")
plot_surface(m_bam_social, "bam_social")
plot_surface(m_bam_rational, "bam_rational")

# ---------------- 8) Variante colapsada a conteos (opcional recomendado) ---
# Colapsa por celda de parámetros (misma especificación estadística)
df_counts_total <- df_all %>%
  transmute(
    IUL, MSP, tau_mean, tau_sd,
    seed_type, contagion_type, network_type,
    successes = ifelse(is.finite(adopt_total), round(adopt_total * 1), NA_real_), # placeholder
    # Vamos a sumar directamente de conteos si están disponibles:
    num_adopters = suppressWarnings(as.numeric(num_adopters)),
    N_nodes      = suppressWarnings(as.numeric(N_nodes_actual))
  )

# Si tenemos conteos originales, los usamos (preferido); si no, salimos
if (all(c("num_adopters","N_nodes_actual") %in% names(df_all))) {
  df_counts_total <- df_all %>%
    transmute(
      IUL, MSP, tau_mean, tau_sd,
      seed_type, contagion_type, network_type,
      successes = as.numeric(num_adopters),
      trials    = as.numeric(N_nodes_actual)
    ) %>%
    group_by(IUL, MSP, tau_mean, tau_sd, seed_type, contagion_type, network_type) %>%
    summarise(
      successes = sum(successes, na.rm = TRUE),
      trials    = sum(trials,    na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    filter(trials > 0)
  
  form_counts <- as.formula(paste(
    "cbind(successes, trials - successes) ~",
    paste(rhs_common, collapse = " + ")
  ))
  
  t0 <- Sys.time()
  m_bam_total_counts <- bam(
    formula  = form_counts,
    family   = binomial,
    data     = df_counts_total,
    method   = "fREML",
    discrete = TRUE,
    nthreads = NUM_CORES_TO_USE
  )
  t1 <- Sys.time()
  message("[BAM Total Counts] elapsed (s): ", round(as.numeric(difftime(t1, t0, units = "secs")), 2))
  
  saveRDS(m_bam_total_counts, file = file.path(OUT_DIR_MODELS, "bam_total_counts.rds"))
  sink(file.path(OUT_DIR_MODELS, "bam_total_counts_summary.txt"))
  on.exit(sink(), add = TRUE)
  print(summary(m_bam_total_counts))
  suppressWarnings(print(gam.check(m_bam_total_counts)))
  png(file.path(OUT_DIR_FIGS, "bam_total_counts_surfaces.png"),
      width = 1600, height = 800, res = 150)
  plot(m_bam_total_counts, scheme = 2, pages = 1, shade = TRUE)
  dev.off()
  
  # Social (si hay conteos y columna)
  df_counts_social <- df_all %>%
    transmute(
      IUL, MSP, tau_mean, tau_sd,
      seed_type, contagion_type, network_type,
      successes = as.numeric(num_adopted_social),
      trials    = as.numeric(N_nodes_actual)
    ) %>%
    group_by(IUL, MSP, tau_mean, tau_sd, seed_type, contagion_type, network_type) %>%
    summarise(
      successes = sum(successes, na.rm = TRUE),
      trials    = sum(trials,    na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    filter(trials > 0)
  
  form_counts_social <- as.formula(paste(
    "cbind(successes, trials - successes) ~",
    paste(rhs_common, collapse = " + ")
  ))
  
  t0 <- Sys.time()
  m_bam_social_counts <- bam(
    formula  = form_counts_social,
    family   = binomial,
    data     = df_counts_social,
    method   = "fREML",
    discrete = TRUE,
    nthreads = NUM_CORES_TO_USE
  )
  t1 <- Sys.time()
  message("[BAM Social Counts] elapsed (s): ", round(as.numeric(difftime(t1, t0, units = "secs")), 2))
  
  saveRDS(m_bam_social_counts, file = file.path(OUT_DIR_MODELS, "bam_social_counts.rds"))
  sink(file.path(OUT_DIR_MODELS, "bam_social_counts_summary.txt"))
  on.exit(sink(), add = TRUE)
  print(summary(m_bam_social_counts))
  suppressWarnings(print(gam.check(m_bam_social_counts)))
  png(file.path(OUT_DIR_FIGS, "bam_social_counts_surfaces.png"),
      width = 1600, height = 800, res = 150)
  plot(m_bam_social_counts, scheme = 2, pages = 1, shade = TRUE)
  dev.off()
  
  # Racional (si hay conteos)
  df_counts_rational <- df_all %>%
    transmute(
      IUL, MSP, tau_mean, tau_sd,
      seed_type, contagion_type, network_type,
      successes = as.numeric(num_adopted_rational),
      trials    = as.numeric(N_nodes_actual)
    ) %>%
    group_by(IUL, MSP, tau_mean, tau_sd, seed_type, contagion_type, network_type) %>%
    summarise(
      successes = sum(successes, na.rm = TRUE),
      trials    = sum(trials,    na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    filter(trials > 0)
  
  form_counts_rational <- as.formula(paste(
    "cbind(successes, trials - successes) ~",
    paste(rhs_common, collapse = " + ")
  ))
  
  t0 <- Sys.time()
  m_bam_rational_counts <- bam(
    formula  = form_counts_rational,
    family   = binomial,
    data     = df_counts_rational,
    method   = "fREML",
    discrete = TRUE,
    nthreads = NUM_CORES_TO_USE
  )
  t1 <- Sys.time()
  message("[BAM Rational Counts] elapsed (s): ", round(as.numeric(difftime(t1, t0, units = "secs")), 2))
  
  saveRDS(m_bam_rational_counts, file = file.path(OUT_DIR_MODELS, "bam_rational_counts.rds"))
  sink(file.path(OUT_DIR_MODELS, "bam_rational_counts_summary.txt"))
  on.exit(sink(), add = TRUE)
  print(summary(m_bam_rational_counts))
  suppressWarnings(print(gam.check(m_bam_rational_counts)))
  png(file.path(OUT_DIR_FIGS, "bam_rational_counts_surfaces.png"),
      width = 1600, height = 800, res = 150)
  plot(m_bam_rational_counts, scheme = 2, pages = 1, shade = TRUE)
  dev.off()
} else {
  message("Conteos originales (num_adopters / N_nodes_actual) no disponibles para colapsar a binomial por celdas. Se omite variante 'counts'.")
}

message("\n== DONE ==")