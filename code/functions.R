# prepare spatial environments
prep_spatial_env <- function() {
  # get map
  uk_mask <- geodata::gadm(country='GBR', level=1, path=".")
  Scottish_border <- uk_mask %>% 
    tidyterra::filter(NAME_1 %in% "Scotland") %>%
    st_as_sf() %>%
    st_transform(crs = "+proj=longlat +datum=WGS84")
  
  Scottish_border_mainland <- Scottish_border %>% 
    st_cast("POLYGON") %>%
    dplyr::mutate(area = st_area(.)) %>%
    arrange(desc(area)) %>%
    slice(1) %>%
    dplyr::select(-area) %>%
    st_transform(crs = 27700)
  
  Scottish_border_mainland <- st_transform(
    Scottish_border_mainland,
    gsub("units=m", "units=km", st_crs(Scottish_border_mainland)$proj4string)
  ) 
  
  # build mesh
  mesh <- fm_mesh_2d(
    boundary = Scottish_border_mainland,
    max.edge = c(0.5, 40),
    cutoff = 10,
    crs = crs(Scottish_border_mainland)
  )
  
  # set up the spde model
  spde <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(150, .5),   # Pr(practic.range<150 km)=0.5
    prior.sigma = c(1, 0.5)     # PR(sigma>1)=0.5
  )
  
  # mesh for prediction
  pxl <- fm_pixels(mesh, mask = Scottish_border_mainland, dims = c(200, 200))
  
  # common objects
  return(list(
    border = Scottish_border_mainland,
    mesh = mesh,
    spde = spde,
    pxl = pxl
  ))
}

# # data preparation: cleaning, rename, factorize
# prep_question_data <- function(df, standard_ans, suffix, border_sf, global_pids) {
#   df_clean <- df %>%
#     dplyr::filter(lat > 0) %>%
#     dplyr::mutate(
#       pid = factor(pid, levels = global_pids), # shared component, keep standard name
#       uni = as.factor(uni),
#       gender = as.factor(gender),
#       age_scaled = as.numeric(scale(age)),
#       region = as.factor(region),
#       bi_ans = ifelse(answer == standard_ans, 0, 1) # 0 for standard, 1 for dialects
#     ) %>%
#     # add suffix
#     dplyr::rename_with(
#       ~ paste0(., "_", suffix), 
#       c(uni, gender, age_scaled, region, bi_ans)
#     ) %>%
#     st_as_sf(coords = c("lng", "lat"), crs = "+proj=longlat +datum=WGS84") %>%
#     st_transform(crs = st_crs(border_sf))
#   
#   # filter points inside mainland Scotland
#   df_mainland <- df_clean[border_sf, ]
#   return(df_mainland)
# }
# data preparation: cleaning, rename, factorize
prep_question_data <- function(df, type, standard_ans = NULL, suffix = NULL, border_sf, global_pids = NULL) {
  # basic factorization
  df_clean <- df %>%
    dplyr::filter(lat > 0) %>%
    dplyr::mutate(
      uni = as.factor(uni),
      gender = as.factor(gender),
      age_scaled = as.numeric(scale(age)),
      region = as.factor(region)
    )
  
  # binarize answer according to question type
  if (type %in% c("sound", "word")) {
    if (is.null(standard_ans)) stop("For 'sound' or 'word' types, standard_ans must be provided.")
    df_clean <- df_clean %>%
      dplyr::mutate(bi_ans = ifelse(answer == standard_ans, 0, 1))     # 0 for standard, 1 for dialects
  } else if (type == "say") {
    df_clean <- df_clean %>%
      dplyr::mutate(bi_ans = ifelse(as.numeric(answer) == 4, 1, 0))    # 1 for usage of dialects, 0 for others (1,2,3)
  } else {
    stop("Type must be 'sound', 'word', or 'say'.")
  }
  
  # pid factorization
  if ("pid" %in% names(df_clean)) {
    if (!is.null(global_pids)) {
      # for joint models
      df_clean <- df_clean %>% dplyr::mutate(pid = factor(pid, levels = global_pids))
    } else {
      # for separate models
      df_clean <- df_clean %>% dplyr::mutate(pid = as.factor(pid))
    }
  }
  
  # add suffix for joint models
  cols_to_rename <- c("uni", "gender", "age_scaled", "region", "bi_ans")
  cols_to_rename <- intersect(cols_to_rename, names(df_clean)) 
  if (!is.null(suffix) && suffix != "") {
    df_clean <- df_clean %>%
      dplyr::rename_with(
        ~ paste0(., "_", suffix),
        dplyr::all_of(cols_to_rename)
      )
  }
  
  # unify coordinate system
  df_sf <- df_clean %>%
    st_as_sf(coords = c("lng", "lat"), crs = "+proj=longlat +datum=WGS84") %>%
    st_transform(crs = st_crs(border_sf))
  
  # truncate on mainland
  df_mainland <- df_sf[border_sf, ]
  
  return(df_mainland)
}


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