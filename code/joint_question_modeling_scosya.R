# This script applies binary classification to 1-5 rating questions 
# and accounts for multiple questions in a joint model using metaprogramming
# set working direction to main folder (".../scottish_dialect")

library(ggplot2)
library(sf)
library(spatstat)
library(inlabru)
library(INLA)
library(SpatialEpi)
library(dplyr)
library(terra)
library(fmesher)
library(scico)
library(patchwork)
library(tidyr)
source("code/functions.R", verbose=FALSE)

# Set category and questions to loop through
category <- "adjectives"  # choose from: "adjectives", "agreement", "auxiliaries", "comp", "have_raising", "imperatives", "inversion", "locative_discovery_expressions", "nominals", "participles", "polarity", "perpositions"    
Q <- c("adjective_shift-e4", "adjective_shift-e5", "intensifiers-q37")   # for "adjectives"
#Q <- c("measures-q11", "measures-q12", "measures-q13", "measures-q14", "measures-q15", "much_many-q5", "much_many-q6", "northern_subject_rule_a12", "northern_subject_rule_a13", "northern_subject_rule_a17", "northern_subject_rule_a29") # for "agreement"   
# other questions ...

# Choose to let the questions share the agegroup covariate, or model separately
share_covariates <- TRUE  

# Calculation of common spatial variables
spatial_env <- prep_spatial_env()
Scottish_border_mainland <- spatial_env$border
mesh_1 <- spatial_env$mesh
spde <- spatial_env$spde
pxl <- spatial_env$pxl

# Data processing
raw_data_list <- list()
all_qids <- c()

# Read, clean, and map rating values
for (i in seq_along(Q)) {
  # Dynamically build filename
  filepath <- sprintf("data/scosya/%s-%s.csv", category, Q[i])

  df <- read.csv(filepath, skip = 1, stringsAsFactors = FALSE) %>% 
    mutate(
      rating = suppressWarnings(as.numeric(rating)),
      lng = suppressWarnings(as.numeric(display_lng)),
      lat = suppressWarnings(as.numeric(display_lat))
    ) %>%
    filter(!is.na(rating) & rating %in% c(1, 2, 4, 5)) %>%
    filter(!is.na(lng) & !is.na(lat)) %>%
    mutate(
      # Binarize rating: 4,5 -> 1 (Dialect); 1,2 -> 0 (Standard); delete 3
      bi_ans = ifelse(rating >= 4, 1, 0),
      agegroup = as.factor(agegroup)
    ) %>%
    distinct(qid, .keep_all = TRUE)  #remove people with multiple answers
  
  raw_data_list[[i]] <- df
  all_qids <- unique(c(all_qids, df$qid))
}

# Convert to SF and handle suffix appending
sf_list <- list()
for (i in seq_along(Q)) {
  suffix <- paste0("q", i)
  df <- raw_data_list[[i]]
  
  df$qid <- factor(df$qid, levels = all_qids)
  df_rename <- df %>%
    rename_with(~ paste0(., "_", suffix), c(bi_ans, agegroup))
  df_sf <- df_rename %>%
    st_as_sf(coords = c("lng", "lat"), crs = "+proj=longlat +datum=WGS84") %>%
    st_transform(crs = st_crs(Scottish_border_mainland))
  
  df_mainland <- df_sf[Scottish_border_mainland, ]
  df_mainland$agegroup <- df_mainland[[paste0("agegroup_", suffix)]]
  sf_list[[i]] <- df_mainland
}
names(sf_list) <- paste0("q", seq_along(Q))


# qid as individual random effect 
cmp_str <- "~ -1 + Intercept_pop(1) + field_pop(geometry, model = spde) + qid_eff(qid, model = 'iid')"

if (share_covariates) {
  cmp_str <- paste(cmp_str, "+ agegroup_eff_shared(agegroup, model = 'iid')")
}

for (i in seq_along(Q)) {
  suffix <- paste0("q", i)
  q_comps <- sprintf("+ Intercept_%s(1) + field_pop_copy_%s(geometry, copy = 'field_pop', fixed = FALSE) + field_%s(geometry, model = spde)", suffix, suffix, suffix)
  
  if (!share_covariates) {
    q_comps <- paste(q_comps, sprintf("+ agegroup_eff_%s(agegroup_%s, model = 'iid')", suffix, suffix))
  }
  
  cmp_str <- paste(cmp_str, q_comps)
}

cmp_multi <- as.formula(cmp_str)


# dynamic likelihood generation
liks_list <- list()

for (i in seq_along(Q)) {
  suffix <- paste0("q", i)
  current_sf <- sf_list[[i]]
  
  # Likelihood for Point Process
  lik_pop <- bru_obs(
    "cp",
    formula = geometry ~ Intercept_pop + field_pop,
    data = current_sf,
    domain = list(geometry = mesh_1),
    samplers = Scottish_border_mainland
  )
  liks_list <- append(liks_list, list(lik_pop))
  
  # Likelihood for Binomial Response
  if (share_covariates) {
    form_str <- sprintf("bi_ans_%s ~ Intercept_%s + agegroup_eff_shared + qid_eff + field_pop_copy_%s + field_%s", suffix, suffix, suffix, suffix)
  } else {
    form_str <- sprintf("bi_ans_%s ~ Intercept_%s + agegroup_eff_%s + qid_eff + field_pop_copy_%s + field_%s", suffix, suffix, suffix, suffix, suffix, suffix)
  }
  
  lik_ans <- bru_obs(
    "binomial",
    formula = as.formula(form_str),
    data = current_sf,
    Ntrials = 1
  )
  liks_list <- append(liks_list, list(lik_ans))
}

# fit the joint model
cat("Fitting full joint model with", length(Q), "questions for category [", category, "]...\n")
t1 <- Sys.time()
fit_multi <- do.call(bru, c(list(components = cmp_multi), liks_list))
t2 <- Sys.time()
cat("Model fitting took", round(difftime(t2, t1, units="secs"), 2), "seconds to complete.\n")


# prediction and visualization
pred_str <- "~ data.frame("
for (i in seq_along(Q)) {
  suffix <- paste0("q", i)
  pred_item <- sprintf("Prob_%s = plogis(Intercept_%s + field_pop_copy_%s + field_%s)", suffix, suffix, suffix, suffix)
  if (i < length(Q)) pred_item <- paste0(pred_item, ", ")
  pred_str <- paste(pred_str, pred_item)
}
pred_str <- paste(pred_str, ")")

pred_response <- predict(fit_multi, pxl, as.formula(pred_str))

plot_list <- list()
for (i in seq_along(Q)) {
  suffix <- paste0("q", i)
  q_name <- Q[i]
  df_prob <- pred_response[[paste0("Prob_", suffix)]]
  
  p_mean <- ggplot(data = df_prob) +
    geom_sf(aes(color = mean), size = 0.5) +  
    scale_color_scico(palette = "roma", direction = -1, limits = c(0, 1), na.value = "grey80") +
    geom_sf(data = Scottish_border_mainland, fill = NA, color = "black", linewidth = 0.2) +
    ggtitle(sprintf("%s (Mean)", q_name)) + labs(color = "Mean") + theme_minimal()
  
  p_sd <- ggplot(data = df_prob) +
    geom_sf(aes(color = sd), size = 0.5) +  
    scale_color_scico(palette = "lajolla", na.value = "grey80") +
    geom_sf(data = Scottish_border_mainland, fill = NA, color = "black", linewidth = 0.2) +
    ggtitle("SD") + labs(color = "SD") + theme_minimal()
  
  plot_list[[i]] <- p_mean + p_sd
}

p_all_combined <- wrap_plots(plot_list, ncol = 1) + 
  plot_annotation(title = sprintf("Joint Model Posterior Predictions: %s (Baseline)", toupper(category)))

# Save the unified plot
out_file <- sprintf("output/figures/scosya_joint_predictions_%s.pdf", category)
ggsave(
  filename = out_file,
  plot = p_all_combined,
  width = 14, height = 4 * length(Q), limitsize = FALSE
)
    