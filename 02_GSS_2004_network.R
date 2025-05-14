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

# n_egos <- nrow(gss_egor$ego) # Necesitarías esto y el peso_mas_pequeno si lo calculas dinámicamente

control_params <- control.ergm.ego(
  # --- Parámetros de la Pseudo-población ---
  ppopsize = 10000, # Tamaño de la pseudo-población para la estimación.
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
    MCMC.burnin = 200000,    # Aumentado para mayor estabilidad con pseudo-poblaciones grandes
    MCMC.interval = 2000,     # Aumentado
    MCMC.samplesize = 2048,  # Número de propuestas por paso.
    seed = 123,              # Para reproducibilidad
    # Podrías necesitar ajustar MPLE.max.dyad.types si tienes muchos atributos categóricos
    # MPLE.max.dyad.types = 1e7 # (ejemplo)
    main.method = "MCMLE"    # Asegurarse de usar MCMLE
  ),
  
  # --- Otros Parámetros ---
  ignore.max.alters = TRUE # (Default) Ignora restricciones en el número máximo de nominaciones de alters.
  # TRUE es el default más reciente y recomendado.
)


# --- Definir popsize (Tamaño de la población de USA en 2004) ---
# Este valor es una aproximación. Debes buscar el valor más preciso.
# Por ejemplo, la población adulta de USA.
US_population_2004 <- 225000000 # Ejemplo, ¡VERIFICAR ESTE NÚMERO!

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

# --- Correr Regresión - Opción A (Educación Categórica Unificada) ---

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
# Siempre revisa la convergencia:
# mcmc.diagnostics(fit_A)

# --- Correr Regresión - Opción B (Educación Numérica Imputada) ---
# Usamos gss_egor_B que tiene 'educ_num_final'
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
# Siempre revisa la convergencia:
# mcmc.diagnostics(fit_B)


# --- Generar una red de 1000 nodos a partir de uno de los modelos ajustados ---
# Usaremos fit_A como ejemplo.
# Para simular una red con atributos de gss_egor, usamos gss_egor como 'basis'.
# 'popsize' aquí es el tamaño de la red a simular.

cat("\n\nSimulando red de 1000 nodos desde fit_A...\n")
# El objeto egor original (gss_egor) se usa como base para los atributos
# de los nodos en la red simulada. Asegúrate de que tenga las variables
# usadas en fit_A (como educ_level_unified).
sim_net_1000 <- simulate.ergm.ego(
  object = fit_A,
  nsim = 1,
  popsize = 1000, # Tamaño de la red a simular
  basis = gss_egor, # Usa los atributos de los egos en gss_egor
  control = control.simulate.ergm.ego(
    SAN = control.san(MCMC.burnin = 10000 * 100, MCMC.interval = 10000), # SAN puede necesitar burn-in largo
    seed = 789
  ),
  verbose = TRUE
)

# sim_net_1000 será una lista de redes (network objects), aquí con 1 elemento
print(sim_net_1000[[1]])

# Puedes inspeccionar la red:
# library(network)
# plot(sim_net_1000[[1]], vertex.cex=0.5)
# summary(sim_net_1000[[1]] ~ edges + nodematch("educ_level_unified")) # etc.