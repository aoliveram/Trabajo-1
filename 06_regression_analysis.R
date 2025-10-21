## ---- Provisional GAM: GRAND_COMBINED-only, robust load, fit --------------

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr); library(stringr); library(mgcv)
})

# ---------------------------------------------------------------------------
# 0) Directories (as you provided)
# ---------------------------------------------------------------------------
dir_atp_emp <- "trabajo_1_files/ATP_diffusion_simulation_files_sigm"
dir_atp_er  <- "trabajo_1_files/ATP_ER_diffusion_simulation_files_sigm"
dir_gss_emp <- "trabajo_1_files/GSS_diffusion_simulation_files_sigm"
dir_gss_er  <- "trabajo_1_files/GSS_ER_diffusion_simulation_files_sigm"

# ---------------------------------------------------------------------------
# 1) List only GRAND_COMBINED .rds files
# ---------------------------------------------------------------------------
ls_grand_combined <- function(path) {
  if (!dir.exists(path)) return(character(0))
  f <- list.files(path, pattern = "\\.rds$", full.names = TRUE, recursive = FALSE)
  f[grepl("GRAND_COMBINED", basename(f), ignore.case = TRUE)]
}

files_atp_emp <- ls_grand_combined(dir_atp_emp)
files_atp_er  <- ls_grand_combined(dir_atp_er)
files_gss_emp <- ls_grand_combined(dir_gss_emp)
files_gss_er  <- ls_grand_combined(dir_gss_er)

# Sanity: stop if nothing found anywhere
if (all(lengths(list(files_atp_emp, files_atp_er, files_gss_emp, files_gss_er)) == 0)) {
  stop("No GRAND_COMBINED .rds files discovered in the four directories.")
}

# ---------------------------------------------------------------------------
# 2) Build the sources table dynamically
# ---------------------------------------------------------------------------
make_sources <- function(paths, contagion_type, network_type) {
  tibble(path = paths,
         contagion_type = contagion_type,
         network_type   = network_type)
}

sources_partial <- bind_rows(
  make_sources(files_atp_emp, "Innovation",       "Empirical"),
  make_sources(files_atp_er,  "Innovation",       "ER"),
  make_sources(files_gss_emp, "CollectiveAction", "Empirical"),
  make_sources(files_gss_er,  "CollectiveAction", "ER")
)

message("Discovered GRAND_COMBINED files:")
print(sources_partial, n = 50)

# ---------------------------------------------------------------------------
# 3) Helpers: parse nesting keys, collect data.frames from nested lists
# ---------------------------------------------------------------------------
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

# Parse seed from filename (works for ATP & GSS names)
# Matches "...means_random.rds" OR "..._random.rds"
parse_seed_from_filename <- function(path) {
  b <- basename(path)
  s1 <- str_match(b, "means_([A-Za-z]+)\\.rds$")[,2]
  s2 <- str_match(b, "_(random|central|eigen|closeness|marginal)\\.rds$")[,2]
  seed <- dplyr::coalesce(s1, s2)
  if (is.na(seed)) seed <- "random"  # fallback
  seed
}

# ---------------------------------------------------------------------------
# 4) Read, label, harmonize one GRAND_COMBINED file
# ---------------------------------------------------------------------------
safe_read_and_label <- function(path, contagion_type, network_type) {
  stopifnot(file.exists(path))
  obj <- readRDS(path)
  
  dfs <- collect_dfs(obj)
  if (!length(dfs)) stop("No data.frames in RDS: ", path)
  
  out <- suppressMessages(bind_rows(dfs))
  
  # Canonical names
  if ("innovation_iul_Gamma" %in% names(out)) out <- rename(out, IUL = innovation_iul_Gamma)
  if ("social_distance_h"    %in% names(out)) out <- rename(out, MSP = social_distance_h)
  
  # Build adoption proportions if counts present (avoid using '.' inside mutate)
  if (!"adopt_total" %in% names(out) && all(c("num_adopters","N_nodes_actual") %in% names(out))) {
    out$adopt_total <- out$num_adopters / out$N_nodes_actual
  }
  if (!"adopt_social" %in% names(out) && all(c("num_adopted_social","N_nodes_actual") %in% names(out))) {
    out$adopt_social <- out$num_adopted_social / out$N_nodes_actual
  }
  if (!"adopt_rational" %in% names(out) && all(c("num_adopted_rational","N_nodes_actual") %in% names(out))) {
    out$adopt_rational <- out$num_adopted_rational / out$N_nodes_actual
  }
  
  # Threshold params if present under other names
  if (!"tau_mean" %in% names(out) && "threshold_mean_param" %in% names(out)) out <- mutate(out, tau_mean = threshold_mean_param)
  if (!"tau_sd"   %in% names(out) && "threshold_sd_param"   %in% names(out)) out <- mutate(out, tau_sd   = threshold_sd_param)
  
  # Seed type (from filename) if missing
  if (!"seed_type" %in% names(out)) out$seed_type <- parse_seed_from_filename(path)
  
  out %>%
    mutate(contagion_type = contagion_type,
           network_type   = network_type)
}

# ---------------------------------------------------------------------------
# 5) Build df_partial (read all discovered GRAND_COMBINED files)
# ---------------------------------------------------------------------------
df_partial <- purrr::pmap_dfr(
  sources_partial,
  function(path, contagion_type, network_type) safe_read_and_label(path, contagion_type, network_type)
)

# Coerce types and drop bad rows
df_partial <- df_partial %>%
  mutate(
    IUL            = as.numeric(IUL),
    MSP            = as.numeric(MSP),
    tau_mean       = as.numeric(tau_mean),
    tau_sd         = as.numeric(tau_sd),
    adopt_total    = as.numeric(adopt_total),
    seed_type      = as.factor(seed_type),
    contagion_type = factor(contagion_type, levels = c("Innovation","CollectiveAction")),
    network_type   = factor(network_type,   levels = c("Empirical","ER"))
  ) %>%
  filter(is.finite(IUL), is.finite(MSP),
         is.finite(tau_mean), is.finite(tau_sd),
         is.finite(adopt_total),
         adopt_total >= 0, adopt_total <= 1) %>%
  tidyr::drop_na(IUL, MSP, tau_mean, tau_sd, seed_type, contagion_type, network_type, adopt_total)

df_partial <- df_partial %>%
  mutate(seed_type = relevel(seed_type, ref = "random"))

# ---------------------------------------------------------------------------
# 6) Coverage report
# ---------------------------------------------------------------------------
message("\n=== Coverage check (GRAND_COMBINED only) ===")
message("N rows: ", nrow(df_partial))
message("\ncontagion_type:\n"); print(table(df_partial$contagion_type, useNA = "ifany"))
message("\nnetwork_type:\n");   print(table(df_partial$network_type,   useNA = "ifany"))
message("\nseed_type:\n");      print(table(df_partial$seed_type,      useNA = "ifany"))

has_multi_seed       <- nlevels(df_partial$seed_type)      >= 2
has_both_contagions  <- nlevels(df_partial$contagion_type) >= 2
has_both_networks    <- nlevels(df_partial$network_type)   >= 2

# ---------------------------------------------------------------------------
# 7) Safe formula (include factors only if they have ≥2 levels)
# ---------------------------------------------------------------------------
smooth_term <- if (has_both_contagions) "s(IUL, MSP, k=10, by = contagion_type)" else "s(IUL, MSP, k=10)"
rhs <- c(
  smooth_term,
  "tau_mean + tau_sd",
  if (has_multi_seed) "seed_type" else NULL,
  if (has_both_networks) "network_type" else NULL
)
base_formula <- as.formula(paste("adopt_total ~", paste(rhs, collapse = " + ")))
message("\nUsing formula:\n", deparse(base_formula))

# ---------------------------------------------------------------------------
# 8) Fit provisional GAM and plot
# ---------------------------------------------------------------------------
ctrl <- gam.control(nthreads = max(1, parallel::detectCores()-1))

time_init <- Sys.time()

m_test <- gam(
  formula = base_formula,
  family  = binomial, method = "REML",
  data    = df_partial,
  control = ctrl
)

time_fin <- Sys.time()
time_total <- difftime(time_fin, time_init, units = "auto")

cat("\n=== Provisional GAM (adopt_total) — summary ===\n")
print(summary(m_test))
cat("\n=== gam.check(m_test) ===\n")
suppressWarnings(gam.check(m_test))

grDevices::dev.new(noRStudioGD = TRUE); plot(m_test, scheme = 2, pages = 1, shade = TRUE)

# Save
if (!dir.exists("models")) dir.create("models", recursive = TRUE)
saveRDS(m_test, "models/m_test_adopt_total.rds")
cat("\nSaved: models/m_test_adopt_total.rds\n")


# ====================== AERS =================================================

NUM_CORES_TO_USE = 8L

time_init <- Sys.time()

m_bam <- bam(
  formula = adopt_total ~ s(IUL, MSP, k=10, by=contagion_type) +
    tau_mean + tau_sd + seed_type + network_type,
  family  = binomial(link="logit"),
  data    = df_partial,
  method  = "fREML",     # faster smoothing
  discrete = TRUE,       # approximate penalties, 10–100× faster for n>1e6
  nthreads = NUM_CORES_TO_USE
)

time_fin <- Sys.time()
time_total <- difftime(time_fin, time_init, units = "auto")

cat("\n=== Provisional GAM (adopt_total) — summary ===\n")
print(summary(m_bam))
cat("\n=== gam.check(m_test) ===\n")
suppressWarnings(gam.check(m_bam))

grDevices::dev.new(noRStudioGD = TRUE); plot(m_test, scheme = 2, pages = 1, shade = TRUE)

# Save
if (!dir.exists("models")) dir.create("models", recursive = TRUE)
saveRDS(m_test, "models/m_test_adopt_total.rds")
cat("\nSaved: models/m_test_adopt_total.rds\n")