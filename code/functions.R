# posterior predictions
plot_and_save_prediction <- function(fit_joint, pxl, border_sf, q, type, exc_threshold,
                                     target_gender = "average", 
                                     target_uni = "average", 
                                     target_age = "average", 
                                     age_mean = NULL, age_sd = NULL) {
  
  # prepare mesh
  pxl_pred <- pxl
  components <- c("Intercept_ans", "field_pop_copy", "field_ans")
  
  if (target_gender != "average") {
    pxl_pred$gender <- as.factor(target_gender)
    components <- c(components, "gender_eff")
  }
  
  if (target_uni != "average") {
    pxl_pred$uni <- as.factor(target_uni)
    components <- c(components, "uni_eff")
  }
  
  if (target_age != "average") {
    if (!is.null(age_mean) && !is.null(age_sd)) {
      pxl_pred$age_scaled <- (as.numeric(target_age) - age_mean) / age_sd
    } else {
      pxl_pred$age_scaled <- as.numeric(target_age)
    }
    components <- c(components, "age_eff")
  }
  
  # dynamic generation
  formula_rhs <- paste(components, collapse = " + ")
  formula_str <- sprintf("~ data.frame(
    Spatial_Field = field_ans, 
    Response = plogis(%s),
    Prob_Exceed = (plogis(%s) > %f)
  )", formula_rhs, formula_rhs, exc_threshold)
  
  pred_formula <- as.formula(formula_str)
  
  posterior <- predict(fit_joint, pxl_pred, pred_formula)
  
  # plot and store figure
  if (target_gender == "average" && target_uni == "average" && target_age == "average") {
    subtitle_text <- "Population Average (Baseline)"
    file_suffix <- "average"
  } else {
    subtitle_text <- sprintf("Gender: %s | Education: %s | Age: %s", 
                             target_gender, target_uni, target_age)
    file_suffix <- sprintf("%s_%s_age%s", target_gender, target_uni, target_age)
  }
  p_mean <- ggplot(data = posterior$Response) +
    geom_sf(aes(color = mean), size = 0.5) +  
    scale_color_scico(palette = "roma" , direction = -1, limits = c(0, 1), na.value = "grey80") +   # fix palette from 0 to 1
    geom_sf(data = border_sf, fill = NA, color = "black", linewidth = 0.2) +
    ggtitle("Dialect Probability (Mean)") +
    labs(color = "Mean") +
    theme_minimal()
  
  p_sd <- ggplot(data = posterior$Response) +
    geom_sf(aes(color = sd), size = 0.5) +  
    scale_color_scico(palette = "lajolla", na.value = "grey80") + 
    geom_sf(data = border_sf, fill = NA, color = "black", linewidth = 0.2) +
    ggtitle("Uncertainty (Standard Deviation)") +
    labs(color = "SD") +
    theme_minimal()
  
  p_map <- p_mean + p_sd + plot_annotation(
    title = sprintf("Question: %s", q),
    subtitle = subtitle_text
  )
  
  filename <- sprintf("output/figures/distribution_and_uncertainty_%s_%s_%s.pdf", type, q, file_suffix)
  ggsave(
    filename = filename,
    plot = p_map,
    width = 12,
    height = 6
  )

  return(posterior)   # for excursion set
}