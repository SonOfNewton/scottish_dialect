# This script applies binary classification to all types of questions and account for multiple questions in joint model
# This script is a test with only two questions
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

# calculation of common variables
spatial_env <- prep_spatial_env()
Scottish_border_mainland <- spatial_env$border
mesh_1 <- spatial_env$mesh
spde <- spatial_env$spde
pxl <- spatial_env$pxl

# test with two questions
q1 <- read.csv("data/csv/sounds-about-right-Q1-down.csv")
q2 <- read.csv("data/csv/sounds-about-right-Q2-more.csv")
# remove duplicate pids within each dataset to ensure each individual contributes only one response per question
q1 <- q1 %>% distinct(pid, .keep_all = TRUE)
q2 <- q2 %>% distinct(pid, .keep_all = TRUE)
# a global set of unique PIDs for consistent factor levels across both datasets
all_pids <- unique(c(q1$pid, q2$pid))

# process data
q1_sf <- prep_question_data(
  df = q1, 
  type = "sound", 
  standard_ans = "noun", 
  suffix = "q1",            
  border_sf = Scottish_border_mainland,
  global_pids = all_pids        
)
q2_sf <- prep_question_data(
  df = q2, 
  type = "sound", 
  standard_ans = "door", 
  suffix = "q2",                
  border_sf = Scottish_border_mainland,
  global_pids = all_pids
)

### model components
# model 1: shared sampling intensity, separate fixed and random effects
cmp_multi <- ~ -1 + 
  # shared field of sampling intensity
  Intercept_pop(1) + 
  field_pop(geometry, model = spde) + 
  # shared individual random effect
  pid_eff(pid, model = "iid") +
  # question 1
  Intercept_q1(1) + 
  uni_eff_q1(uni_q1, model = "iid") +          
  gender_eff_q1(gender_q1, model = "iid") +    
  age_eff_q1(age_scaled_q1, model = "linear") +    
  field_pop_copy_q1(geometry, copy = "field_pop", fixed = FALSE) +    #W0
  field_q1(geometry, model = spde) +   #W1
  # question 2
  Intercept_q2(1) + 
  uni_eff_q2(uni_q2, model = "iid") +          
  gender_eff_q2(gender_q2, model = "iid") +    
  age_eff_q2(age_scaled_q2, model = "linear") +    
  field_pop_copy_q2(geometry, copy = "field_pop", fixed = FALSE) +    #W0
  field_q2(geometry, model = spde)     #W1



### likelihoods
# sampling intensity of Q1 and Q2 are considered as two realizations of the same point process
# likelihood 1: Q1 locations
lik_pop_q1 <- bru_obs(
  "cp",
  formula = geometry ~ Intercept_pop + field_pop,
  data = q1_sf,
  domain = list(geometry = mesh_1),
  samplers = Scottish_border_mainland
)
# likelihood 2: Q2 locations
lik_pop_q2 <- bru_obs(
  "cp",
  formula = geometry ~ Intercept_pop + field_pop,
  data = q2_sf,
  domain = list(geometry = mesh_1),
  samplers = Scottish_border_mainland
)

# likelihood 3: binomial response for Q1 
lik_q1 <- bru_obs(
  "binomial",
  formula = bi_ans_q1 ~ Intercept_q1 + uni_eff_q1 + gender_eff_q1 + age_eff_q1 + pid_eff + field_pop_copy_q1 + field_q1,
  data = q1_sf,
  Ntrials = 1
)
# likelihood 4: binomial response for Q2
lik_q2 <- bru_obs(
  "binomial",
  formula = bi_ans_q2 ~ Intercept_q2 + uni_eff_q2 + gender_eff_q2 + age_eff_q2 + pid_eff + field_pop_copy_q2 + field_q2,
  data = q2_sf,
  Ntrials = 1
)


### fit the joint model
fit_multi <- bru(cmp_multi, lik_pop_q1, lik_pop_q2, lik_q1, lik_q2)
# summary(fit_multi)


### predictions and visualization
# plot latent fields (W0, W1, W2)
pred_fields <- predict(
  fit_multi, 
  pxl, 
  ~ data.frame(
    W0_Sampling = field_pop, 
    W1_Q1_Specific = field_q1,
    W2_Q2_Specific = field_q2
  )
)

p_w0 <- ggplot(data = pred_fields$W0_Sampling) +
  geom_sf(aes(color = mean), size = 0.5) +  
  scale_color_scico(palette = "lajolla", na.value = "grey80") +
  ggtitle("W0: Shared Sampling Intensity") + theme_minimal()

p_w1 <- ggplot(data = pred_fields$W1_Q1_Specific) +
  geom_sf(aes(color = mean), size = 0.5) +  
  scale_color_scico(palette = "broc", na.value = "grey80") +
  ggtitle("W1: Q1 (down) Specific Field") + theme_minimal()

p_w2 <- ggplot(data = pred_fields$W2_Q2_Specific) +
  geom_sf(aes(color = mean), size = 0.5) +  
  scale_color_scico(palette = "broc", na.value = "grey80") +
  ggtitle("W2: Q2 (more) Specific Field") + theme_minimal()

p_fields_combined <- p_w0 | p_w1 | p_w2
ggsave(
  filename = "output/figures/joint_spatial_fields_Q1_Q2.pdf",
  plot = p_fields_combined,
  width = 15, height = 10
)


# plot posterior predictions (mean, std) of a baseline profile (omits pid_eff, uni_eff, age_eff, and gender_eff)
pred_response <- predict(
  fit_multi, 
  pxl, 
  ~ data.frame(
    Q1_Prob = plogis(Intercept_q1 + field_pop_copy_q1 + field_q1),
    Q2_Prob = plogis(Intercept_q2 + field_pop_copy_q2 + field_q2)
  )
)

p_q1_mean <- ggplot(data = pred_response$Q1_Prob) +
  geom_sf(aes(color = mean), size = 0.5) +  
  scale_color_scico(palette = "roma", direction = -1, limits = c(0, 1), na.value = "grey80") +
  geom_sf(data = Scottish_border_mainland, fill = NA, color = "black", linewidth = 0.2) +
  ggtitle("Q1 (down) Probability (Mean)") + labs(color = "Mean") + theme_minimal()

p_q1_sd <- ggplot(data = pred_response$Q1_Prob) +
  geom_sf(aes(color = sd), size = 0.5) +  
  scale_color_scico(palette = "lajolla", na.value = "grey80") +
  geom_sf(data = Scottish_border_mainland, fill = NA, color = "black", linewidth = 0.2) +
  ggtitle("Q1 Uncertainty (SD)") + labs(color = "SD") + theme_minimal()

p_q2_mean <- ggplot(data = pred_response$Q2_Prob) +
  geom_sf(aes(color = mean), size = 0.5) +  
  scale_color_scico(palette = "roma", direction = -1, limits = c(0, 1), na.value = "grey80") +
  geom_sf(data = Scottish_border_mainland, fill = NA, color = "black", linewidth = 0.2) +
  ggtitle("Q2 (more) Probability (Mean)") + labs(color = "Mean") + theme_minimal()

p_q2_sd <- ggplot(data = pred_response$Q2_Prob) +
  geom_sf(aes(color = sd), size = 0.5) +  
  scale_color_scico(palette = "lajolla", na.value = "grey80") +
  geom_sf(data = Scottish_border_mainland, fill = NA, color = "black", linewidth = 0.2) +
  ggtitle("Q2 Uncertainty (SD)") + labs(color = "SD") + theme_minimal()

p_response_combined <- (p_q1_mean + p_q1_sd) / (p_q2_mean + p_q2_sd) + 
  plot_annotation(title = "Joint Model Posterior Predictions (Population Average Baseline)")

ggsave(
  filename = "output/figures/joint_predictions_Q1_Q2.pdf",
  plot = p_response_combined,
  width = 14, height = 12
)

