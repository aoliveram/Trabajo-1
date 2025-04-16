# --- Notes ---

# The sript collapses a lot of info about 'race' (eg. 'Asian' --> 'Other'), 'relig' (eg. 'Christian' --> 'Other'), etc.

# --- 0. Setup ---
# Install packages if you haven't already
# install.packages("tidyverse")
# install.packages("egor")
# install.packages("remotes")
# remotes::install_github("statnet/ergm.ego")
# Load libraries
#library(tidyverse)
library(tidyr)
library(dplyr)
library(egor)
library(ergm)
library(ergm.ego)

# Set seed for reproducibility
set.seed(12345)

# --- 1. Load Data ---
gss_wide_data_raw <- read.csv("trabajo_1_files/GSS_2004_EGO_bin.csv")

# --- 2. Prepare Ego Data ---
# (Ego data preparation remains the same as the previous version)
# Add a unique ego ID
gss_wide_data <- gss_wide_data_raw %>%
  mutate(.egoID = row_number()) %>%
  select(.egoID, everything())

# Define ego attributes
ego_cols <- c(".egoID", "sex", "race", "educ", "age", "relig", "degree", "numgiven")

# Create the egos data frame, applying codebook info
egos_df <- gss_wide_data %>%
  select(all_of(ego_cols)) %>%
  mutate(
    across(c(educ, relig), ~na_if(., 99)),
    degree = na_if(degree, 8)
  ) %>%
  mutate(
    sex_ego = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    race_ego = factor(race, levels = c(1, 2, 3), labels = c("White", "Black", "Other")),
    degree_ego = factor(degree,
                        levels = 0:4,
                        labels = c("Less than HS", "High school", "Associate/Junior college", "Bachelor's", "Graduate")),
    relig_ego_recoded = factor(case_when(
      relig == 1 ~ "Protestant",
      relig == 2 ~ "Catholic",
      relig == 3 ~ "Jewish",
      relig == 4 ~ "None",
      relig %in% c(5:13) ~ "Other",
      TRUE ~ NA_character_
    ), levels=c("Protestant", "Catholic", "Jewish", "None", "Other")),
    race_ego_recoded = factor(case_when(
      race == 1 ~ "White",
      race == 2 ~ "Black",
      race == 3 ~ "Other",
      TRUE ~ NA_character_
    ), levels=c("White", "Black", "Other"))
  ) %>%
  filter(!is.na(numgiven) & numgiven != 9) %>%
  select(.egoID, numgiven, sex_ego, race_ego, race_ego_recoded, educ, age, relig_ego_recoded, degree_ego)


# --- 3. Prepare Alter Data (Handling Status/Talkto NAs) ---

# Define stems
alter_attr_stems <- c("sex", "race", "educ", "age", "relig")
alter_status_stems <- c("spouse", "parent", "sibling", "child", "othfam", "cowork", "memgrp", "neighbr", "friend", "advisor", "other")
alter_talkto_stems <- c("talkto")
stems_to_pivot <- c(alter_attr_stems, alter_status_stems, alter_talkto_stems)
cols_to_pivot <- intersect(unlist(lapply(stems_to_pivot, function(stem) paste0(stem, 1:5))), names(gss_wide_data))

# Stems where 9 = NA (based on codebook for alters)
stems_9_is_na <- c("sex", "race", "educ", "relig", # From previous check
                   "spouse", "parent", "sibling", "child", "othfam", # Status vars
                   "cowork", "memgrp", "neighbr", "friend", "advisor", "other", # Status vars
                   "talkto") # Talkto var

# Stems where 8 = NA/DK
stems_8_is_na <- c("educ", "relig") # From previous check

# Generate column names for NA handling
cols_9_is_na <- intersect(unlist(lapply(stems_9_is_na, function(s) paste0(s, 1:5))), names(gss_wide_data))
cols_8_is_na <- intersect(unlist(lapply(stems_8_is_na, function(s) paste0(s, 1:5))), names(gss_wide_data))

# Prepare a temporary wide df for pivoting, handling NAs *before* pivoting
gss_wide_for_pivot <- gss_wide_data %>%
  mutate(across(all_of(cols_9_is_na), ~na_if(., 9)), # Convert 9 to NA for relevant vars
         across(all_of(cols_8_is_na), ~na_if(., 8))  # Convert 8 to NA for relevant vars
  )

# Pivot longer - capturing original variable name and value
alters_long <- gss_wide_for_pivot %>%
  select(.egoID, numgiven, all_of(cols_to_pivot)) %>%
  pivot_longer(
    cols = all_of(cols_to_pivot),
    names_to = "original_var",
    values_to = "value",
    values_drop_na = FALSE
  ) %>%
  extract(original_var, into = c("stem", ".altID"), regex = "([a-zA-Z]+)([1-5])", remove = FALSE, convert = TRUE) %>%
  filter(.altID <= numgiven) %>%
  select(.egoID, .altID, stem, value)

# Recode values *before* pivoting wider
alters_long_recoded <- alters_long %>%
  mutate(
    value_char = as.character(value),
    recoded_value = case_when(
      stem == "race" & value == 4 ~ "White",
      stem == "race" & value == 2 ~ "Black",
      stem == "race" & value %in% c(1, 3, 5) ~ "Other",
      stem == "relig" & value == 1 ~ "Protestant",
      stem == "relig" & value == 2 ~ "Catholic",
      stem == "relig" & value == 3 ~ "Jewish",
      stem == "relig" & value == 4 ~ "None",
      stem == "relig" & value == 5 ~ "Other",
      TRUE ~ value_char
    ),
    stem = case_when(
      stem == "race" ~ "race_alt_recoded",
      stem == "relig" ~ "relig_alt_recoded",
      TRUE ~ stem
    )
  ) %>%
  select(.egoID, .altID, stem, recoded_value)

# Pivot wider using the 'stem' column for names and 'recoded_value' for values
alters_wide_again <- alters_long_recoded %>%
  group_by(.egoID, .altID, stem) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_wider(
    names_from = stem,
    values_from = recoded_value
  )

# --- Process Statuses and Finalize Alters DF (with talkto as factor) ---
status_vars_in_wide <- intersect(alter_status_stems, names(alters_wide_again))

alters_df <- alters_wide_again %>%
  mutate(across(all_of(status_vars_in_wide), ~as.numeric(.))) %>%
  mutate(relationship = factor(case_when(
    spouse == 1 ~ "Spouse",
    parent == 1 ~ "Parent",
    sibling == 1 ~ "Sibling",
    child == 1 ~ "Child",
    othfam == 1 ~ "Other Family",
    cowork == 1 ~ "Coworker",
    memgrp == 1 ~ "Group Member",
    neighbr == 1 ~ "Neighbor",
    friend == 1 ~ "Friend",
    advisor == 1 ~ "Advisor",
    other == 1 ~ "Other",
    TRUE ~ "Unknown"
  ))) %>%
  mutate(
    sex = as.numeric(sex),
    educ = as.numeric(educ),
    age = as.numeric(age),
    talkto = as.numeric(talkto), # Convert talkto to numeric first
    
    sex_alt = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    race_alt_recoded = factor(race_alt_recoded, levels=c("White", "Black", "Other")),
    relig_alt_recoded = factor(relig_alt_recoded, levels=c("Protestant", "Catholic", "Jewish", "None", "Other")),
    educ_alt = factor(educ,
                      levels = 0:7,
                      labels = c("1-6 years", "7-9 years", "10-12 years", "HS grad", "Some college", "Associate", "Bachelor's", "Graduate")),
    # Create ORDERED factor for talkto
    talkto_alt = factor(talkto,
                        levels = 1:4, # Based on codebook (1-4)
                        labels = c("Every day", "Weekly", "Monthly", "Less often"),
                        ordered = TRUE) # Specify it's ordered
  ) %>%
  # Select final columns, including the new talkto factor
  select(
    .egoID, .altID,
    sex_alt, race_alt_recoded, educ_alt, age, relig_alt_recoded,
    relationship, talkto_alt # Use the factor talkto_alt
  )

# --- 4. Filter Data based on NAs in Modeling Variables ---
# Add talkto_alt to the list if using it in the model
model_vars_ego <- c("age", "sex_ego", "race_ego_recoded", "degree_ego", "relig_ego_recoded")
model_vars_alter <- c("age", "sex_alt", "race_alt_recoded", "educ_alt", "relig_alt_recoded", "relationship", "talkto_alt") # Added talkto_alt

# (Filtering code remains the same)
alters_df_filtered <- alters_df %>%
  filter(complete.cases(select(., all_of(model_vars_alter))))
valid_egos_post_filter <- unique(alters_df_filtered$.egoID)

egos_df_filtered <- egos_df %>%
  filter(.egoID %in% valid_egos_post_filter) %>%
  filter(complete.cases(select(., all_of(model_vars_ego))))
valid_egos_post_filter <- intersect(valid_egos_post_filter, egos_df_filtered$.egoID)

egos_df_final <- egos_df_filtered %>% filter(.egoID %in% valid_egos_post_filter)
alters_df_final <- alters_df_filtered %>% filter(.egoID %in% valid_egos_post_filter)
gss_wide_data_final <- gss_wide_data %>% filter(.egoID %in% valid_egos_post_filter)

print(paste("Number of egos after final filtering:", nrow(egos_df_final)))
print(paste("Number of alters after final filtering:", nrow(alters_df_final)))
print("Final Alter Data Frame Head:")
print(head(alters_df_final))
print("Final Alter Data Frame Structure:")
str(alters_df_final)
print("Alter Talkto Frequencies:")
print(table(alters_df_final$talkto_alt, useNA = "ifany"))


# --- 5. Create egor Object ---
# (Code remains the same)

print(head(egos_df_final))
print(head(alters_df_final))

# ERRORS START -- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

gss_egor <- egor(
  egos = egos_df_final,
  alters = alters_df_final,
  aaties = gss_wide_data_final,
  ID.vars = list(ego = ".egoID", alter = ".altID", source = ".egoID")#,
  #aaties.vars = list(name.vars = "close", name.sep = "", ego.id = ".egoID")
)

# --- 6. Examine the egor Object ---
# (Code remains the same)
print("egor object summary:")
summary(gss_egor)

# --- 7. Define and Fit Ego-ERGM ---
# Add nodefactor("talkto_alt") to the model
model_formula <- S ~ edges +
  ego("age") +
  ego("sex_ego") +
  ego("degree_ego") +
  nodeattr("age") +
  nodefactor("sex_alt") +
  nodefactor("educ_alt") +
  nodefactor("relationship") +
  nodefactor("talkto_alt") + # Added talkto factor effect
  absdiff("age") +
  nodematch("sex_ego", "sex_alt") +
  nodematch("race_ego_recoded", "race_alt_recoded") +
  nodematch("relig_ego_recoded", "relig_alt_recoded")

# (Fitting code remains the same)
print("Fitting Ego-ERGM model...")
control_settings <- control.ergm.ego(MCMLE.maxit = 40, seed = 123)
fit <- tryCatch({
  ergm.ego(model_formula, object = gss_egor, control = control_settings)
}, error = function(e) {
  print(paste("Error during model fitting:", e$message))
  NULL
})

# --- 8 & 9. Simulate Networks & Examine ---
# (Code remains the same as previous version)
if (!is.null(fit)) {
  print("Model Summary:")
  summary(fit)
  nsim <- 1000
  print(paste("Simulating", nsim, "networks..."))
  simulated_networks <- simulate(fit, nsim = nsim, control = control.simulate.ergm.ego(seed = 456))
  print(paste("Generated", length(simulated_networks), "simulated networks."))
  if (length(simulated_networks) > 0) {
    print("Density of first 5 simulated networks:")
    print(sapply(simulated_networks[1:min(5, nsim)], network::network.density))
  }
} else {
  print("Skipping simulation due to model fitting error.")
}

print("Script finished.")