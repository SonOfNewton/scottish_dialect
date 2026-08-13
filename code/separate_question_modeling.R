# This script applies binary classification to all types of questions and models each question separately
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

# "word": "give your word" maps 
# "sound": "sounds about right" maps 
# "say": "I would never say that" maps 
type = "word"  
use_factors = TRUE   # incoporate age, gender, and education factor
plot_average = FALSE   # plot posterior mean with an average age, gender, and education
# only used when plot_average = FALSE 
specify_age = 20
specify_gender = "Female"
specify_uni = "Y"

plot_fields = FALSE   # can choose to skip the plot of spatial field w1 and w2
do_excursion = FALSE   # can choose to skip the excursion set plot
exc_threshold = 0.6        # threshold of dialect usage proportion (each type should have different threshold)
exc_confidence = 0.85      # threshold of confidence

# questions
if (type == "sound"){
  Q <- c("Q1-down", "Q2-more","Q4-die", "Q5-head", "Q6-make","Q7-soft", "Q8-home", "Q9-card", "Q10-stone")
  standard <- c("noun", "door", "tie", "bed", "bake", "croft", "foam", "scarred", "loan")
} else if (type =="say"){
  Q <- c("Q1-Are-you-wanting-to-come-with-me", "Q2-The-cat-needs-fed", "Q3-You-telt-me-that-already","Q4-Im-coming-with-you,-amnt-I", 
         "Q8-Im-going-to-my-bed", "Q20-What-are-youse-doing-tonight", "Q26-When-are-you-back-at-the-school")
} else if (type == "word"){
  Q <- c("Q1-knock", "Q8-fed-up", "Q11-remember", "Q13-yes", "Q16-cry", "Q18-church", "Q23-known", "Q28-halloween")
  standard <- c("knock", "fed up", "remember", "Yes", "crying", "church", "known", "trick-or-treating")
}

all_uni_summary <- data.frame()    # for education effect
all_age_summary <- data.frame()    # for age effect
all_gender_summary <- data.frame() # for gender effect
hyper_table <- data.frame()        # for posterior of hyperparameters
spatial_fields_df <- data.frame()  # for spatial field

# calculation of common variables
spatial_env <- prep_spatial_env()
Scottish_border_mainland <- spatial_env$border
mesh_1 <- spatial_env$mesh
spde <- spatial_env$spde
pxl <- spatial_env$pxl

# model components
if (use_factors == TRUE){
  cmp_joint <- ~ -1 + 
    # items in spatial point process
    Intercept_pop(1) + 
    field_pop(geometry, model = spde) + 
    # items of independent variables
    Intercept_ans(1) + 
    uni_eff(uni, model = "iid") +          # intercept of education
    gender_eff(gender, model = "iid") +    # intercept of gender
    age_eff(age_scaled, model = "linear") +    # age effect (scaled)
    #region_eff(region, model = "iid") +    # intercept of region
    field_pop_copy(geometry, copy = "field_pop", fixed = FALSE) + # shared spatial field (times \alpha)
    field_ans(geometry, model = spde)       # residual for answers
} else{
  cmp_joint <- ~ -1 + 
    # items in spatial point process
    Intercept_pop(1) + 
    field_pop(geometry, model = spde) + 
    # items of independent variables
    Intercept_ans(1) + 
    field_pop_copy(geometry, copy = "field_pop", fixed = FALSE) + # shared spatial field (times \alpha)
    field_ans(geometry, model = spde)       # residual for answers
}


### iterate over all questions
for (i in seq_along(Q)) {
  q <- Q[i]
  if (type == "sound"){
    s <- standard[i]
    question <- read.csv(sprintf("data/sfy/sounds-about-right-%s.csv", q))
  } else if(type == "say"){
    question <- read.csv(sprintf("data/sfy/grammar-%s.csv", q))
  } else if (type == "word"){
    s <- standard[i]
    question <- read.csv(sprintf("data/sfy/lexical-%s.csv", q))
  }
  question <- question[question$lat > 0, ]

  # process data
  sf_mainland_answers <- prep_question_data(
    df = question, 
    type = type, 
    standard_ans = s, 
    suffix = NULL,       
    global_pids = NULL,  
    border_sf = Scottish_border_mainland
  )
  
  sf_mainland <- sf_mainland_answers %>% dplyr::select(geometry)

  #likelihoods
  lik_pop <- bru_obs(
    "cp",
    formula = geometry ~ Intercept_pop + field_pop,
    data = sf_mainland,
    domain = list(geometry = mesh_1),
    samplers = Scottish_border_mainland
  )
  
  if (use_factors == TRUE){
    lik_ans <- bru_obs(
      "binomial",
      formula = bi_ans ~ Intercept_ans + uni_eff + gender_eff + age_eff + field_pop_copy + field_ans,
      data = sf_mainland_answers,
      Ntrials = 1
    )
  } else{
    lik_ans <- bru_obs(
      "binomial",
      formula = bi_ans ~ Intercept_ans + field_pop_copy + field_ans,
      data = sf_mainland_answers,
      Ntrials = 1
    )
  }
  
  fit_joint <- bru(cmp_joint, lik_pop, lik_ans)
  
  # extract hyperparameters
  hyper <- fit_joint$summary.hyperpar
  hyper_out <- hyper %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Question = q) %>%
    relocate(Question, Parameter)
  
  hyper_table <- bind_rows(hyper_table, hyper_out)
  
  # extract education effect
  if ("uni_eff" %in% names(fit_joint$summary.random)) {
    uni_summary <- fit_joint$summary.random$uni_eff
    uni_summary$Question <- q
    all_uni_summary <- bind_rows(all_uni_summary, uni_summary)
  }
  
  # extract gender effect
  if ("gender_eff" %in% names(fit_joint$summary.random)) {
    gender_summary <- fit_joint$summary.random$gender_eff
    gender_summary$Question <- q
    all_gender_summary <- bind_rows(all_gender_summary, gender_summary)
  }
  
  # extract age effect
  if ("age_eff" %in% rownames(fit_joint$summary.fixed)) {
    age_summary <- fit_joint$summary.fixed["age_eff", , drop=FALSE] 
    age_summary <- tibble::rownames_to_column(age_summary, "Parameter")
    age_summary$Question <- q
    all_age_summary <- bind_rows(all_age_summary, age_summary)
  }
  
  ### visualization
  # reverse logit
  posterior <- predict(
    fit_joint, 
    pxl, 
    ~ data.frame(
      Spatial_Field = field_ans, 
      Response = plogis(Intercept_ans + field_pop_copy + field_ans),
      Prob_Exceed = (plogis(Intercept_ans + field_pop_copy + field_ans) > exc_threshold)  # for excursion set
    )
  )
  

  if (plot_average == TRUE){
    # for an average profile
    posterior <- plot_and_save_prediction(
      fit_joint = fit_joint, 
      pxl = pxl, 
      border_sf = Scottish_border_mainland, 
      q = q, 
      type = type, 
      exc_threshold = exc_threshold
    ) 
  } else {
    current_age_mean <- mean(answers_full$age, na.rm = TRUE)
    current_age_sd <- sd(answers_full$age, na.rm = TRUE)
    # for a specific profile (e.g.: a 50-year-old man that did not go to university)
    plot_and_save_prediction(
      fit_joint = fit_joint, 
      pxl = pxl, 
      border_sf = Scottish_border_mainland, 
      q = q, 
      type = type, 
      exc_threshold = exc_threshold,
      target_gender = specify_gender,    
      target_uni = specify_uni,      
      target_age = specify_age,           
      age_mean = current_age_mean, 
      age_sd = current_age_sd
    )
  }
  
  
  # posterior mean of spatial random field
  pred_spatial <- predict(
    fit_joint, 
    pxl, 
    ~ data.frame(Pop_Density_Field = field_pop, Answer_Specific_Field = field_ans)
  )
  
  # store spatial field
  if (nrow(spatial_fields_df) == 0) {
    spatial_fields_df <- data.frame(row_id = 1:nrow(pred_spatial$Answer_Specific_Field))
  }
  spatial_fields_df[[q]] <- pred_spatial$Answer_Specific_Field$mean
  
  if (plot_fields==TRUE){
    # for W_pop
    p2 <- ggplot(data = pred_spatial$Pop_Density_Field) +
      geom_sf(aes(color = mean), size = 0.5) +  
      scale_color_scico(palette = "lajolla", na.value = "grey80") +
      ggtitle("Latent field for point process") +
      theme_minimal()
    
    # for W_ans
    p3 <- ggplot(data = pred_spatial$Answer_Specific_Field) +
      geom_sf(aes(color = mean), size = 0.5) +  
      scale_color_scico(palette = "broc", na.value = "grey80") +
      ggtitle("Independent latent field for answer") +
      theme_minimal()
    
    p_spatial <- p2 + p3
    ggsave(
      filename = sprintf("output/figures/spatial_fields_%s_%s.pdf", type, q),
      plot = p_spatial,
      width = 12,
      height = 6
    ) 
  }
  
  ### excursion plot
  if (do_excursion == TRUE){
    # Left panel: P( p(s) > 0.6 | data )
    p_exc_prob <- ggplot(data = posterior$Prob_Exceed) +
      geom_sf(aes(color = mean), size = 0.5) +  
      scale_color_scico(palette = "oslo", direction = -1, limits = c(0, 1), na.value = "grey80") +
      geom_sf(data = Scottish_border_mainland, fill = NA, color = "black", linewidth = 0.2) +
      ggtitle(sprintf("P( p(s) > %g | data )", exc_threshold)) +
      labs(color = "Probability") +
      theme_minimal()
    
    # Right panel: Excursion Set (Red region where P > 0.85)
    posterior$Prob_Exceed <- posterior$Prob_Exceed %>%
      mutate(is_confident = ifelse(mean > exc_confidence, "Yes", "No"))
    p_exc_set <- ggplot(data = posterior$Prob_Exceed) +
      geom_sf(aes(color = is_confident), size = 0.5) +  
      scale_color_manual(values = c("Yes" = "#d73027", "No" = "grey80")) +
      geom_sf(data = Scottish_border_mainland, fill = NA, color = "black", linewidth = 0.2) +
      ggtitle(sprintf("Confidence > %g%% for p(s) > %g", exc_confidence*100, exc_threshold)) +
      labs(color = "Excursion Set") +
      theme_minimal()
    
    p_excursion_combined <- p_exc_prob + p_exc_set + plot_annotation(title = sprintf("Marginal Excursion Analysis: %s", q))
    ggsave(
      filename = sprintf("output/figures/excursion_set_%s_%s.pdf", type, q),
      plot = p_excursion_combined,
      width = 12,
      height = 6
    )
  }
}


### analysis
# education effect plot
if(nrow(all_uni_summary) > 0) {
  p_uni <- ggplot(
    all_uni_summary,
    aes(
      x = Question,
      y = mean,
      ymin = `0.025quant`,
      ymax = `0.975quant`,
      color = ID
    )
  ) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +  
    labs(
      x = "Question",
      y = "Effect on log-odds",
      color = "Education"
    )
  ggsave(
    filename = sprintf("output/figures/education_effect_%s.pdf", type),
    plot = p_uni,
    width = 12,
    height = 6
  )
}

# gender effect plot
if(nrow(all_gender_summary) > 0) {
  p_gender <- ggplot(
    all_gender_summary,
    aes(
      x = Question,
      y = mean,
      ymin = `0.025quant`,
      ymax = `0.975quant`,
      color = ID
    )
  ) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +  
    labs(
      x = "Question",
      y = "Effect on log-odds",
      color = "Gender"
    )
  ggsave(
    filename = sprintf("output/figures/gender_effect_%s.pdf", type),
    plot = p_gender,
    width = 12,
    height = 6
  )
}

# age effect plot
if(nrow(all_age_summary) > 0) {
  p_age <- ggplot(
    all_age_summary,
    aes(
      x = Question,
      y = mean,
      ymin = `0.025quant`,
      ymax = `0.975quant`
    )
  ) +
    geom_pointrange(color = "#d73027", size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +  
    labs(
      title = sprintf("Effect of Age on Dialect Probability [%s]", type),
      subtitle = "Coefficient for 1 Standard Deviation increase in Age",
      x = "Question",
      y = "Effect on log-odds (Coefficient)"
    )
  ggsave(
    filename = sprintf("output/figures/age_effect_%s.pdf", type),
    plot = p_age,
    width = 10,
    height = 6
  )
}

tablename <- sprintf("output/tables/hyperparameters_%s.csv", type)
write.csv(hyper_table, tablename, row.names = FALSE)  

tablename <- sprintf("output/tables/spatial_field_%s.csv", type)
write.csv(spatial_fields_df, tablename, row.names = FALSE)  

# plot for alpha
hyper_table <- read.csv(sprintf("output/tables/hyperparameters_%s.csv", type), check.names = FALSE)
alpha_summary <- hyper_table %>%
  filter(grepl("field_pop_copy", Parameter))

p_alpha <- ggplot(
  alpha_summary,
  aes(
    x = Question,
    y = mean,
    ymin = `0.025quant`,
    ymax = `0.975quant`
  )
) +
  geom_pointrange(color = "black", size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = sprintf("Shared Spatial Field Scaling Factor (Alpha) by Question [%s]", type),
    subtitle = "Posterior Mean and 95% Credible Intervals",
    x = "Question",
    y = "Alpha (Effect size)"
  )

ggsave(
  filename = sprintf("output/figures/alpha_%s.pdf", type),
  plot = p_alpha,
  width = 10,
  height = 6
)

# similarity analysis of questions
spatial_fields_df <- read.csv(sprintf("output/tables/spatial_field_%s.csv", type))
spatial_data_only <- spatial_fields_df %>% dplyr::select(-row_id)

# Pearson correlation
cor_matrix <- cor(spatial_data_only, use = "pairwise.complete.obs", method = "pearson")
write.csv(cor_matrix, sprintf("output/tables/spatial_correlation_matrix_%s.csv", type))

dist_matrix <- as.dist(1 - cor_matrix)    # range: from 0 (most similar) to 2 (most different)
hc <- hclust(dist_matrix, method = "ward.D2")    # hierarchical clustering

# plot
pdf(sprintf("output/figures/dendrogram_spatial_clustering_%s.pdf", type), width = 12, height = 8)
plot(hc, 
     main = sprintf("Hierarchical Clustering of Dialect Words (%s)", type),
     ylab = "Distance (1 - Pearson r)",
     cex = 0.9,     
     hang = -1)    
dev.off()


