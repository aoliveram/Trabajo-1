networks_dir <- "trabajo_1_files/ATP_network_ergm/"
graphs_ATP <- list()

for (i in 1:2) {
  ATP_net <- readRDS(paste0(networks_dir, "ATP_network_simulated_1000_mur_", sprintf("%03d", i), ".rds"))
  ATP_net <- asIgraph(ATP_net)
  
  graphs_ATP[[i]] <- ATP_net
}

current_graph_obj_sim <- graphs_ATP[[1]]

# Calcular la matriz de disimilitud de Gower
attributes_for_distance <- data.frame(
  age = V(current_graph_obj_sim)$age,
  educ_num = V(current_graph_obj_sim)$educ_num,
  race = as.factor(V(current_graph_obj_sim)$race), # de 'character' a 'factor'
  relig = as.factor(V(current_graph_obj_sim)$relig),
  sex = as.factor(V(current_graph_obj_sim)$sex)
)
social_distance_matrix_d_ij <- as.matrix(daisy(attributes_for_distance, metric = "gower"))

hist(social_distance_matrix_d_ij[lower.tri(social_distance_matrix_d_ij)])

# ¿Cómo afecta cada variable al puntaje total?

n <- nrow(attributes_for_distance)

# Crear todas las combinaciones de pares únicos (i < j)
pairs_df <- as.data.frame(t(combn(n, 2)))
colnames(pairs_df) <- c("i", "j")

# Añadir la distancia Gower para cada par
pairs_df$d_ij <- social_distance_matrix_d_ij[lower.tri(social_distance_matrix_d_ij)]

# Agregar las diferencias por variable
pairs_df <- pairs_df %>%
  mutate(
    age_diff = abs(attributes_for_distance$age[i] - attributes_for_distance$age[j]),
    educ_diff = abs(attributes_for_distance$educ_num[i] - attributes_for_distance$educ_num[j]),
    race_diff = as.integer(attributes_for_distance$race[i] != attributes_for_distance$race[j]),
    relig_diff = as.integer(attributes_for_distance$relig[i] != attributes_for_distance$relig[j]),
    sex_diff = as.integer(attributes_for_distance$sex[i] != attributes_for_distance$sex[j])
  )

# Corremos regresión
model <- lm(d_ij ~ age_diff + educ_diff + race_diff + relig_diff + sex_diff, data = pairs_df)
summary(model)

# Cada variable categórica aporta con ~0.2 a la distancia social.
# Coefficients:
#   Estimate Std. Error   t value Pr(>|t|)    
# (Intercept) 4.270e-13  1.827e-15 2.338e+02   <2e-16 ***
#   age_diff    2.817e-03  4.419e-17 6.374e+13   <2e-16 ***
#   educ_diff   1.000e-02  2.778e-16 3.600e+13   <2e-16 ***
#   race_diff   2.000e-01  1.284e-15 1.558e+14   <2e-16 ***
#   relig_diff  2.000e-01  1.422e-15 1.407e+14   <2e-16 ***
#   sex_diff    2.000e-01  1.260e-15 1.588e+14   <2e-16 ***