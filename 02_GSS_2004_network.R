# Crearemos la red. Los datos de 
# gss_egor <- readRDS("trabajo_1_files/gss_egor.rds")
# gss_egos <- gss_egor$ego
# gss_alters <- gss_egor$alter
# pueden ser cargados, pero para EDUC se debe correr 
#
# 01_GSS_2004_alter_educ.R
# 
# Eso genera la columna "educ_level_unified". 

library(ergm.ego)
library(dplyr)

# --- Cálculo de ppopsize dinámico ---
n_egos_real <- nrow(gss_egor$ego)

# Normalizar pesos para que la media sea aproximadamente 1
# Esto es importante para la interpretación de ppopsize en relación con n_egos
# y para evitar problemas numéricos si los pesos son muy grandes o muy pequeños.
# Si ya están normalizados, esto no debería cambiar mucho.
mean_weight <- mean(gss_egor$ego$wtssnr, na.rm = TRUE)
gss_egor$ego$wtssnr_normalized <- gss_egor$ego$wtssnr / mean_weight

# Calcular el peso mínimo (excluyendo NAs o ceros si los hubiera)
min_normalized_weight <- min(gss_egor$ego$wtssnr_normalized[gss_egor$ego$wtssnr_normalized > 0], na.rm = TRUE)

# Krivitsky & Morris (2017) sugieren al menos 1 observación en la pseudo-población
# por cada "tipo" de ego después de ponderar. Multiplicar por 3 es un factor de seguridad.
ppopsize_calculated <- ceiling(max(1000, (n_egos_real / min_normalized_weight) * 1)) # Usamos *1 según K&M

cat("Tamaño de Pseudo-población Calculado:", ppopsize_calculated, "\n")

# --- Definir popsize (Tamaño de la población de USA en 2004) ---

US_population_2004 <- 225000000 # Tentativo, verificar

# --- Parámetros de Control ---

# NO usado
init_coeffs_temporal <- c(edges = -19.0, # Valor de ejemplo
                          `nodematch.race.BLACK` = 1.5, # Nombre y valor de ejemplo
                          `nodematch.race.OTHER` = 1.5, # Nombre y valor de ejemplo
                          `nodematch.sex.FEMALE` = 0.3,  # Nombre y valor de ejemplo
                          absdiff.age = -0.05,       # Valor de ejemplo
                          `gwesp.fixed.0.25` = 0.1)   # Valor inicial para gwesp

control_params <- control.ergm.ego(
  # --- Parámetros de la Pseudo-población ---
  ppopsize = ppopsize_calculated, # Tamaño de la pseudo-población para la estimación.
  # Debería ser suficientemente grande, especialmente con datos ponderados.
  # Un valor más preciso sería: max(1000, n_egos / min_weight * 3)
  # ppopsize.mul es 1 por defecto, así que ppopsize es el valor absoluto.
  
  ppop.wt = "round", # Método para asignar pesos fraccionales de egos en la pseudo-población.
  # "round" (default) redondea al entero más cercano. "sample" remuestrea.
  # "round" es generalmente bueno.
  
  # --- Parámetros de los Estadísticos Suficientes ---
  stats.wt = "data", # Cómo se ponderan las contribuciones de los egos a los estadísticos.
  # "data" (default) usa los pesos originales del diseño muestral.
  # "ppop" usa los pesos finales de la pseudo-población. "data" es más estándar.
  
  stats.est = "asymptotic", # Método para estimar la varianza de los estadísticos.
  # "asymptotic" (default, método Delta) es común.
  # Otras opciones son "bootstrap", "jackknife", "survey", "naive".
  
  # --- Parámetros para la llamada interna a ergm() ---
  ergm = control.ergm(
    # Parámetros de MCMC para el ajuste del modelo ERGM en la pseudo-población
    # Estos son valores iniciales razonables, pero SIEMPRE revisa la convergencia
    # con mcmc.diagnostics(). ergm.ego puede usar defaults más grandes.
    MCMC.burnin = 100000,    # Reducido inicialmente para pruebas más rápidas
    MCMC.interval = 1000,      # Reducido inicialmente
    MCMC.samplesize = 1024,   # Mantenido o ligeramente reducido
    seed = 123,
    # Podrías necesitar ajustar MPLE.max.dyad.types si tienes muchos atributos categóricos
    # MPLE.max.dyad.types = 1e7 # (ejemplo)
    main.method = "MCMLE",    # Asegurarse de usar MCMLE
    #init = init_coeffs_temporal, # <--- init VA AQUÍ, DENTRO DE control.ergm()
    MCMLE.maxit = 60 # Mantener o ajustar si es necesario
    # eval.loglik = FALSE # Puedes probar esto si sigue habiendo problemas iniciales
  ),
  
  # --- Otros Parámetros ---
  ignore.max.alters = TRUE, # (Default) Ignora restricciones en el número máximo de nominaciones de alters.
  # TRUE es el default más reciente y recomendado.
)

# --- ESTRUCTURA DE LA FÓRMULA PARA ergm.ego ---
#
# 1. DENSIDAD BASE:
#    edges: Controla la propensión general a formar lazos.
#
# 2. HOMOFILIA POR EDUCACIÓN:
#    Opción A (Categórica - usando 'educ_level_unified'):
#      nodematch("educ_level_unified"): Tendencia a lazos con mismo nivel educativo.
#      nodefactor("educ_level_unified"):  Diferencias en actividad por nivel.
#    Opción B (Numérica - usando 'educ_num_final'):
#      absdiff("educ_num_final"): Probabilidad de lazo disminuye con diferencia en años de educ.
#
# 3. HOMOFILIA POR RAZA (usando 'race'):
#    nodematch("race"): Tendencia a lazos con misma raza.
#    nodefactor("race"):  Diferencias en actividad por raza.
#
# 4. HOMOFILIA POR SEXO (usando 'sex'):
#    nodematch("sex"): Tendencia a lazos con mismo sexo.
#    nodefactor("sex"):  Diferencias en actividad por sexo.
#
# 5. HOMOFILIA POR EDAD (usando 'age'):
#    absdiff("age"): Probabilidad de lazo disminuye con diferencia de edad.
#    nodecov("age"):  Efecto lineal de la edad en la actividad.
#
# 6. HOMOFILIA POR RELIGIÓN (usando 'relig'):
#    nodematch("relig"): Tendencia a lazos con misma religión.
#    nodefactor("relig"):  Diferencias en actividad por religión.
#
# 7. ESTRUCTURA DE RED (TRANSITIVIDAD/AGRUPAMIENTO):
#    gwesp(alpha, fixed=TRUE): Geometrically Weighted Edgewise Shared Partners.
#                               Captura la tendencia de "amigos de amigos son amigos".
#                               alpha es un parámetro de decaimiento. 0.5 es un valor inicial.
#                               Aumentar alpha da más peso a configuraciones más densas.
#
# (Otros términos opcionales podrían ser degree(0) para aislados, etc.)
#

# --- Opción A (Educación Categórica Unificada) ---

# --- Versión básica ---
fit_a <- ergm.ego(gss_egor ~ edges +
                    nodematch("educ_level_unified") +
                    nodematch("race") +
                    nodematch("sex") +
                    absdiff("age") +
                    nodematch("relig"),
                  popsize = US_population_2004,
                  control = control_params,
                  verbose = TRUE)

summary(fit_a)
mcmc.diagnostics(fit_a)

# --- Versión con GWESP ---
fit_A <- ergm.ego(gss_egor ~ edges +
                    nodematch("educ_level_unified") +
                    nodematch("race") +
                    nodematch("sex") +
                    absdiff("age") +
                    nodematch("relig") +
                    gwesp(0.5, fixed = TRUE),
                  popsize = US_population_2004,
                  control = control_params,
                  verbose = TRUE) # O FALSE para menos output

summary(fit_A)
# mcmc.diagnostics(fit_A)

# --- Opción B (Educación Numérica Imputada) ---
cat("\n\nAjustando modelo Opción B (educación numérica)...\n")
fit_B <- ergm.ego(gss_egor_B ~ edges +
                    absdiff("educ_num_final") + # Usando la variable numérica
                    nodematch("race") +
                    nodematch("sex") +
                    absdiff("age") +
                    nodematch("relig") +
                    gwesp(0.5, fixed = TRUE),
                  popsize = US_population_2004,
                  control = control_params,
                  verbose = TRUE) # O FALSE para menos output

summary(fit_B)
# mcmc.diagnostics(fit_B)


# --- Generar una red de 1000 nodos a partir de uno de los modelos ajustados ---

# El objeto egor original (gss_egor) se usa como base para los atributos
# de los nodos en la red simulada. Asegúrate de que tenga las variables
# usadas en fit_A (como educ_level_unified).
sim_net_1000 <- simulate(
  fit_a,
  popsize = 1000, # Tamaño de la red a simular
  nsim = 1,
  control = control.simulate.ergm.ego(
    simulate=control.simulate.formula(MCMC.burnin=2e6)
  ),
  verbose = TRUE
)


print(sim_net_1000[[1]])

library(network)
plot(sim_net_1000[[1]], vertex.cex=0.5)
summary(sim_net_1000[[1]] ~ edges + nodematch("educ_level_unified")) # etc.
