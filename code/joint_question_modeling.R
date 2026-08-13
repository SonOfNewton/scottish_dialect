# This script applies binary classification to all types of questions and account for multiple questions in joint model
# This script applies metaprogramming techniques to build a model that contains an arbitrary number of terms (questions)
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
source("code/functions.R",verbose=FALSE)

share_covariates <- TRUE   # choose to let the questions share covariates for age, gender, and education, or to model separately for each question

# calculation of common variables
spatial_env <- prep_spatial_env()
Scottish_border_mainland <- spatial_env$border
mesh_1 <- spatial_env$mesh
spde <- spatial_env$spde
pxl <- spatial_env$pxl

# all questions
Q <- c("Q1-down", "Q2-more", "Q4-die", "Q5-head")#, "Q6-make", "Q7-soft", "Q8-home", "Q9-card", "Q10-stone")
standard <- c("noun", "door", "tie", "bed")#, "bake", "croft", "foam", "scarred", "loan")
# note that computational cost for model fitting grows exponentially with number of question included

raw_data_list <- list()
all_pids <- c()

for (i in seq_along(Q)) {
  df <- read.csv(sprintf("data/sfy/sounds-about-right-%s.csv", Q[i])) %>% 
    distinct(pid, .keep_all = TRUE)
  raw_data_list[[i]] <- df
  all_pids <- unique(c(all_pids, df$pid))
}

sf_list <- list()
for (i in seq_along(Q)) {
  suffix <- paste0("q", i)
  # process data
  df_sf <- prep_question_data(
    df = raw_data_list[[i]], 
    type = "sound", 
    standard_ans = standard[i], 
    suffix = suffix,            
    border_sf = Scottish_border_mainland,
    global_pids = all_pids        
  )
  
  # for model 1
  df_sf$uni <- df_sf[[paste0("uni_", suffix)]]
  df_sf$gender <- df_sf[[paste0("gender_", suffix)]]
  df_sf$age_scaled <- df_sf[[paste0("age_scaled_", suffix)]]
  
  sf_list[[i]] <- df_sf
}
names(sf_list) <- paste0("q", seq_along(Q))


# basic shared components
cmp_str <- "~ -1 + Intercept_pop(1) + field_pop(geometry, model = spde) + pid_eff(pid, model = 'iid')"

if (share_covariates) {
  # shared fixed and random effects
  cmp_str <- paste(cmp_str, "+ uni_eff_shared(uni, model = 'iid') + gender_eff_shared(gender, model = 'iid') + age_eff_shared(age_scaled, model = 'linear')")
}

# question-specific components
for (i in seq_along(Q)) {
  suffix <- paste0("q", i)
  q_comps <- sprintf("+ Intercept_%s(1) + field_pop_copy_%s(geometry, copy = 'field_pop', fixed = FALSE) + field_%s(geometry, model = spde)", suffix, suffix, suffix)
  
  if (!share_covariates) {
    q_comps <- paste(q_comps, sprintf("+ uni_eff_%s(uni_%s, model = 'iid') + gender_eff_%s(gender_%s, model = 'iid') + age_eff_%s(age_scaled_%s, model = 'linear')", suffix, suffix, suffix, suffix, suffix, suffix))
  }
  
  cmp_str <- paste(cmp_str, q_comps)
}

# transform string to formula
cmp_multi <- as.formula(cmp_str)

### likelihoods
liks_list <- list()

for (i in seq_along(Q)) {
  suffix <- paste0("q", i)
  current_sf <- sf_list[[i]]
  
  # point process for sampling intensity
  lik_pop <- bru_obs(
    "cp",
    formula = geometry ~ Intercept_pop + field_pop,
    data = current_sf,
    domain = list(geometry = mesh_1),
    samplers = Scottish_border_mainland
  )
  liks_list <- append(liks_list, list(lik_pop))
  
  # point process for binary response
  if (share_covariates) {
    form_str <- sprintf("bi_ans_%s ~ Intercept_%s + uni_eff_shared + gender_eff_shared + age_eff_shared + pid_eff + field_pop_copy_%s + field_%s", suffix, suffix, suffix, suffix)
  } else {
    form_str <- sprintf("bi_ans_%s ~ Intercept_%s + uni_eff_%s + gender_eff_%s + age_eff_%s + pid_eff + field_pop_copy_%s + field_%s", suffix, suffix, suffix, suffix, suffix, suffix, suffix)
  }
  
  lik_ans <- bru_obs(
    "binomial",
    formula = as.formula(form_str),
    data = current_sf,
    Ntrials = 1
  )
  liks_list <- append(liks_list, list(lik_ans))
}


### fit the joint model
cat("Fitting full joint model with", length(Q), "questions...\n")
t1 <- Sys.time()
# equivalent to: bru(cmp_multi, liks_list[[1]], liks_list[[2]], ..., liks_list[[18]])
fit_multi <- do.call(bru, c(list(components = cmp_multi), liks_list))
t2 <- Sys.time()
cat("Model fitting took", t2-t1, "seconds to complete.")
# summary(fit_multi)

### predictions and visualization
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
  plot_annotation(title = "joint model posterior predictions (all questions baseline)")

ggsave(
  filename = "output/figures/joint_predictions_all.pdf",
  plot = p_all_combined,
  width = 14, height = 4 * length(Q), limitsize = FALSE
)
