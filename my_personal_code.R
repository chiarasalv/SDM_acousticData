# ============================================================
# 02_run_PLS_revision_analysis_TOT.R
# REVISION PIPELINE — Picus viridis SDM
# PLS ONLY — RESAMPLED 90 m + RAO
# MAIN ANALYSIS = TOT
# ============================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(terra)
library(caret)
library(pROC)
library(pls)
library(ggplot2)
library(purrr)
library(viridis)
library(ape)
library(tibble)

set.seed(123)


# 1) PATHS AND INPUTS
setwd("my/path")

prepared_dir <- "prepared_resampled90_data"
out_dir <- "outputs_revision_PLS_TOT_resampled90_rao"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

main_period <- "TOT"

df <- read.csv(file.path(
  prepared_dir,
  "BIRDS_SPECIES_ALLMETRICS_RESAMPLED90_WITH_RAO.csv"
))

ahr_metrics  <- rast(file.path(prepared_dir, "ahr_metrics_resampled90_masked.tif"))
slud_metrics <- rast(file.path(prepared_dir, "slud_metrics_resampled90_masked.tif"))
vals_metrics <- rast(file.path(prepared_dir, "vals_metrics_resampled90_masked.tif"))

logger_points <- read.csv(file.path(
  prepared_dir,
  "logger_coordinates.csv"
))

area_rasters <- list(
  ahr  = ahr_metrics,
  slud = slud_metrics,
  vals = vals_metrics
)


# 2) SETTINGS
target_species <- "Picus viridis"

pr_start <- as.Date("2024-03-01")
pr_end   <- as.Date("2024-06-20")

env_vars <- c(
  "vegetation_cover", "canopy_volume", "chm_rugosity", "num_trees",
  "zmax", "zmean", "zsd", "zskew", "zkurt",
  "pzabovezmean", "pzabove2",
  paste0("zq", seq(5, 95, 5)),
  paste0("zpcum", 1:9),
  "NDVI_mean", "NDVI_sd", "NDVI_max",
  "rao_chm", "rao_ndvi"
)

env_vars_map <- env_vars


# 3) OPTIONAL: LOGGER REMOVAL
# vals_to_remove <- c("vals_moth9", "vals_moth4")
# 
# df <- df %>%
#   filter(!id_area %in% vals_to_remove)
# 
# logger_points <- logger_points %>%
#   filter(!id_area %in% vals_to_remove)


# 4) CLEAN DATA
df_clean <- df %>%
  mutate(
    Date = as.Date(Date),
    biotope = sub("_.*", "", id_area),
    logger = sub(".*_", "", id_area),
    period = case_when(
      Date >= pr_start & Date <= pr_end ~ "PR",
      TRUE ~ "NR"
    )
  )

print(table(df_clean$biotope))
print(table(df_clean$scientific_name == target_species))


# 5) LOGGER × DAY DATASET
daily_data <- df_clean %>%
  group_by(id_area, biotope, logger, Date, period) %>%
  summarise(
    presence = as.integer(any(scientific_name == target_species)),
    n_detections_total = n(),
    n_species_detected = n_distinct(scientific_name),
    mean_confidence = mean(confidence, na.rm = TRUE),
    .groups = "drop"
  )

env_logger <- df_clean %>%
  select(id_area, all_of(env_vars)) %>%
  distinct()

env_check <- env_logger %>%
  count(id_area) %>%
  filter(n > 1)

print(env_check)

model_data <- daily_data %>%
  left_join(env_logger, by = "id_area") %>%
  drop_na()

data_TOT <- model_data
data_PR  <- model_data %>% filter(period == "PR")
data_NR  <- model_data %>% filter(period == "NR")

data_main <- data_TOT

print(table(model_data$presence))
print(table(model_data$period, model_data$presence))

logger_summary <- model_data %>%
  group_by(biotope, id_area) %>%
  summarise(
    n_days = n(),
    presences = sum(presence),
    absences = sum(presence == 0),
    prevalence = mean(presence),
    mean_detections_total = mean(n_detections_total, na.rm = TRUE),
    mean_species_detected = mean(n_species_detected, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  logger_summary,
  file.path(out_dir, "logger_summary_presence_effort.csv"),
  row.names = FALSE
)


# 6) SITE-LEVEL DATASET — TOT
site_data <- data_main %>%
  group_by(id_area, biotope) %>%
  summarise(
    presence_rate = mean(presence),
    n_days = n(),
    mean_detections_total = mean(n_detections_total, na.rm = TRUE),
    mean_species_detected = mean(n_species_detected, na.rm = TRUE),
    across(all_of(env_vars), first),
    .groups = "drop"
  ) %>%
  left_join(
    logger_points %>% select(id_area, x, y),
    by = "id_area"
  )




# CREATE DATASET LOGGER-LEVEL (ONLY VALS)
vals_logger <- model_data %>%
  filter(biotope == "vals") %>%
  group_by(id_area, biotope) %>%
  summarise(
    n_days = n(),
    presences = sum(presence),
    absences = sum(presence == 0),
    prevalence = mean(presence),
    mean_detections_total = mean(n_detections_total, na.rm = TRUE),
    mean_species_detected = mean(n_species_detected, na.rm = TRUE),
    across(all_of(env_vars), first),
    .groups = "drop"
  ) %>%
  left_join(
    logger_points %>% select(id_area, x, y),
    by = "id_area"
  )

vals_logger


# 2) DISTANZE SPAZIALI TRA LOGGER VALS
coords <- vals_logger %>%
  select(x, y) %>%
  as.matrix()

spatial_dist <- as.matrix(dist(coords))
rownames(spatial_dist) <- vals_logger$id_area
colnames(spatial_dist) <- vals_logger$id_area

spatial_dist


# 3) DISTANZE AMBIENTALI TRA LOGGER VALS
# prediTTORI standardizzati
env_scaled <- vals_logger %>%
  select(all_of(env_vars)) %>%
  scale()


env_dist <- as.matrix(dist(env_scaled))
rownames(env_dist) <- vals_logger$id_area
colnames(env_dist) <- vals_logger$id_area

env_dist


# 4) TABELLA PAIRWISE: DISTANZA SPAZIALE + AMBIENTALE
pairs <- expand.grid(
  logger_1 = vals_logger$id_area,
  logger_2 = vals_logger$id_area
) %>%
  filter(logger_1 < logger_2)

vals_pairwise_redundancy <- pairs %>%
  rowwise() %>%
  mutate(
    spatial_distance_m = round(spatial_dist[logger_1, logger_2], 1),
    environmental_distance = round(env_dist[logger_1, logger_2], 3)
  ) %>%
  ungroup() %>%
  arrange(spatial_distance_m, environmental_distance)

vals_pairwise_redundancy


# 5) PER OGNI LOGGER: VICINO SPAZIALE PIÙ PROSSIMO E LOGGER SIMILI
diag(spatial_dist) <- NA
diag(env_dist) <- NA

vals_redundancy_summary <- lapply(vals_logger$id_area, function(id) {
  
  nearest_spatial <- names(which.min(spatial_dist[id, ]))
  nearest_env <- names(which.min(env_dist[id, ]))
  
  data.frame(
    id_area = id,
    prevalence = vals_logger$prevalence[vals_logger$id_area == id],
    presences = vals_logger$presences[vals_logger$id_area == id],
    absences = vals_logger$absences[vals_logger$id_area == id],
    mean_detections_total = vals_logger$mean_detections_total[vals_logger$id_area == id],
    mean_species_detected = vals_logger$mean_species_detected[vals_logger$id_area == id],
    
    nearest_spatial_logger = nearest_spatial,
    spatial_distance_m = round(spatial_dist[id, nearest_spatial], 1),
    
    nearest_environmental_logger = nearest_env,
    environmental_distance = round(env_dist[id, nearest_env], 3)
  )
}) %>%
  bind_rows() %>%
  arrange(spatial_distance_m, environmental_distance)

vals_redundancy_summary


# 6) decision table
vals_decision_table <- vals_redundancy_summary %>%
  mutate(
    decision = case_when(
      id_area %in% c("vals_moth4", "vals_moth9") ~ "removed",
      TRUE ~ "retained"
    ),
    reason = case_when(
      id_area == "vals_moth4" ~ "spatially/environmentally redundant within Valsura",
      id_area == "vals_moth9" ~ "spatially/environmentally redundant within Valsura",
      TRUE ~ "retained for spatial/environmental representation"
    )
  )

vals_decision_table


# gruppi di ridondanza
vals_decision_table <- vals_decision_table %>%
  mutate(
    redundancy_group = case_when(
      id_area %in% c("vals_moth3", "vals_moth4") ~ "G1",
      id_area %in% c("vals_moth8", "vals_moth9") ~ "G2",
      TRUE ~ NA_character_
    )
  )

vals_decision_table_final <- vals_decision_table %>%
  select(
    redundancy_group,
    id_area,
    decision,
    prevalence,
    mean_detections_total,
    mean_species_detected,
    nearest_spatial_logger,
    spatial_distance_m,
    nearest_environmental_logger,
    environmental_distance,
    reason
  ) %>%
  arrange(
    redundancy_group,
    decision
  )

vals_decision_table_final

colnames(vals_decision_table_final) <- c(
  "Redundancy group",
  "Logger",
  "Decision",
  "Prevalence",
  "Mean detections/day",
  "Mean species/day",
  "Nearest spatial logger",
  "Spatial distance (m)",
  "Nearest environmental logger",
  "Environmental distance",
  "Reason"
)

write.csv2(
  vals_decision_table_final,
  file.path(out_dir,
            "Table_S1_valsura_logger_redundancy_assessment.csv"),
  row.names = FALSE
)

vals_to_remove <- c("vals_moth9", "vals_moth4")

df <- df %>%
  filter(!id_area %in% vals_to_remove)

logger_points <- logger_points %>%
  filter(!id_area %in% vals_to_remove)


# 7) CORRELATION FILTER
filter_predictors <- function(data, predictors, cutoff = 0.75) {
  
  x <- data %>%
    select(all_of(predictors)) %>%
    select(where(is.numeric))
  
  x <- x[, sapply(x, function(v) sd(v, na.rm = TRUE) > 0), drop = FALSE]
  
  if (ncol(x) < 2) {
    return(list(
      kept = names(x),
      removed = character(0),
      cor_mat = NULL
    ))
  }
  
  cor_mat <- cor(x, use = "pairwise.complete.obs")
  remove_idx <- caret::findCorrelation(cor_mat, cutoff = cutoff)
  
  kept <- names(x)[-remove_idx]
  removed <- names(x)[remove_idx]
  
  list(
    kept = kept,
    removed = removed,
    cor_mat = cor_mat
  )
}


# 8) SITE-LEVEL LOOCV WITH PLS — TOT
run_site_loocv_pls <- function(site_data, predictors) {
  
  dat <- site_data %>%
    select(id_area, biotope, presence_rate, x, y, all_of(predictors)) %>%
    drop_na()
  
  preds <- rep(NA, nrow(dat))
  
  for (i in seq_len(nrow(dat))) {
    
    train <- dat[-i, ]
    test  <- dat[i, ]
    
    ncomp_i <- min(3, length(predictors), nrow(train) - 1)
    if (ncomp_i < 1) next
    
    ctrl <- trainControl(method = "none")
    
    set.seed(123)
    
    mod <- train(
      presence_rate ~ .,
      data = train %>% select(-id_area, -biotope, -x, -y),
      method = "pls",
      trControl = ctrl,
      tuneGrid = data.frame(ncomp = ncomp_i)
    )
    
    preds[i] <- predict(mod, newdata = test)
  }
  
  out <- dat %>%
    mutate(
      predicted = preds,
      residual = presence_rate - predicted
    )
  
  list(
    predictions = out,
    cor = cor(out$presence_rate, out$predicted, use = "complete.obs"),
    rmse = sqrt(mean((out$presence_rate - out$predicted)^2, na.rm = TRUE))
  )
}

pred_site <- filter_predictors(site_data, env_vars_map, cutoff = 0.75)

site_loocv_pls <- run_site_loocv_pls(
  site_data = site_data,
  predictors = pred_site$kept
)

write.csv(
  site_loocv_pls$predictions,
  file.path(out_dir, "site_level_LOOCV_predictions_PLS_TOT.csv"),
  row.names = FALSE
)

p_site_loocv <- ggplot(
  site_loocv_pls$predictions,
  aes(x = presence_rate, y = predicted, color = biotope)
) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    x = "Observed full-year acoustic occurrence rate",
    y = "Predicted occurrence rate",
    color = "Biotope",
    title = "Site-level leave-one-logger-out validation — PLS, TOT"
  ) + theme(
    axis.text.x = element_text(size = 12),  # valori asse X
    axis.text.y = element_text(size = 12),  # valori asse Y
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )

ggsave(
  file.path(out_dir, "Figure_1_site_level_LOOCV_PLS_TOT.png"),
  p_site_loocv,
  width = 7,
  height = 5,
  dpi = 300
)


# 9) MORAN'S I ON SITE-LEVEL RESIDUALS
run_moran_residuals <- function(pred_df) {
  
  dat <- pred_df %>%
    filter(!is.na(x), !is.na(y), !is.na(residual))
  
  coords <- as.matrix(dat[, c("x", "y")])
  d <- as.matrix(dist(coords))
  diag(d) <- NA
  
  weights <- 1 / d
  weights[is.na(weights)] <- 0
  weights[is.infinite(weights)] <- 0
  
  moran <- ape::Moran.I(dat$residual, weight = weights)
  
  tibble(
    observed = moran$observed,
    expected = moran$expected,
    sd = moran$sd,
    p_value = moran$p.value
  )
}

moran_site <- run_moran_residuals(site_loocv_pls$predictions)

write.csv(
  moran_site,
  file.path(out_dir, "moran_residuals_site_level_LOOCV_PLS_TOT.csv"),
  row.names = FALSE
)


# 10) LEAVE-ONE-AREA-OUT PLS — TOT
run_leave_one_area_pls <- function(site_data, predictors) {
  
  dat <- site_data %>%
    select(id_area, biotope, presence_rate, all_of(predictors)) %>%
    drop_na()
  
  areas <- unique(dat$biotope)
  all_preds <- data.frame()
  
  for (area in areas) {
    
    train <- dat %>% filter(biotope != area)
    test  <- dat %>% filter(biotope == area)
    
    ncomp_i <- min(3, length(predictors), nrow(train) - 1)
    if (ncomp_i < 1) next
    
    ctrl <- trainControl(method = "none")
    
    mod <- train(
      presence_rate ~ .,
      data = train %>% select(-id_area, -biotope),
      method = "pls",
      trControl = ctrl,
      tuneGrid = data.frame(ncomp = ncomp_i)
    )
    
    pred <- predict(mod, newdata = test)
    
    tmp <- test %>%
      mutate(
        predicted = pred,
        residual = presence_rate - predicted,
        left_out_area = area
      )
    
    all_preds <- bind_rows(all_preds, tmp)
  }
  
  list(
    predictions = all_preds,
    cor = cor(all_preds$presence_rate, all_preds$predicted, use = "complete.obs"),
    rmse = sqrt(mean((all_preds$presence_rate - all_preds$predicted)^2, na.rm = TRUE))
  )
}

area_transfer_pls <- run_leave_one_area_pls(
  site_data = site_data,
  predictors = pred_site$kept
)

write.csv(
  area_transfer_pls$predictions,
  file.path(out_dir, "leave_one_area_out_predictions_PLS_TOT.csv"),
  row.names = FALSE
)

p_area_transfer <- ggplot(
  area_transfer_pls$predictions,
  aes(x = presence_rate, y = predicted, color = left_out_area)
) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    x = "Observed full-year acoustic occurrence rate",
    y = "Predicted occurrence rate",
    color = "Left-out area",
    title = "Leave-one-area-out transferability test — PLS, TOT"
  ) + theme(
    axis.text.x = element_text(size = 12),  # valori asse X
    axis.text.y = element_text(size = 12),  # valori asse Y
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )

ggsave(
  file.path(out_dir, "Figure_2_leave_one_area_out_PLS_TOT.png"),
  p_area_transfer,
  width = 7,
  height = 5,
  dpi = 300
)


# 11) DAILY PLS SDM FUNCTION
run_daily_pls_sdm <- function(data, predictors) {
  
  dat <- data %>%
    select(id_area, biotope, presence, all_of(predictors)) %>%
    drop_na()
  
  dat$presence <- factor(
    dat$presence,
    levels = c(0, 1),
    labels = c("absence", "presence")
  )
  
  if (length(unique(dat$presence)) < 2) {
    message("Skipping PLS: only one response class.")
    return(NULL)
  }
  
  folds <- groupKFold(
    group = dat$id_area,
    k = length(unique(dat$id_area))
  )
  
  ctrl <- trainControl(
    method = "cv",
    index = folds,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    sampling = NULL
  )
  
  tune_grid <- expand.grid(
    ncomp = 1:min(5, length(predictors))
  )
  
  res <- tryCatch({
    
    set.seed(123)
    
    mod <- train(
      presence ~ .,
      data = dat %>% select(-id_area, -biotope),
      method = "pls",
      metric = "ROC",
      trControl = ctrl,
      tuneGrid = tune_grid
    )
    
    pred <- mod$pred
    
    if (!"presence" %in% names(pred)) {
      message("Skipping PLS: probability column missing.")
      return(NULL)
    }
    
    pred <- pred %>% filter(!is.na(presence))
    
    if (length(unique(pred$obs)) < 2) {
      message("Skipping PLS: only one class in predictions.")
      return(NULL)
    }
    
    roc_obj <- pROC::roc(
      response = pred$obs,
      predictor = pred$presence,
      levels = c("absence", "presence")
    )
    
    auc_value <- as.numeric(pROC::auc(roc_obj))
    
    best_thr <- pROC::coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = "threshold"
    )
    
    list(
      model = mod,
      predictions = pred,
      auc = auc_value,
      youden_threshold = as.numeric(best_thr),
      predictors = predictors,
      data = dat
    )
    
  }, error = function(e) {
    message("PLS failed: ", e$message)
    return(NULL)
  })
  
  return(res)
}


# 12) PERIOD COMPARISON — DIAGNOSTIC ONLY
period_datasets <- list(
  TOT = data_TOT,
  PR  = data_PR,
  NR  = data_NR
)

all_results <- list()

for (period_name in names(period_datasets)) {
  
  dat_period <- period_datasets[[period_name]]
  
  filtered <- filter_predictors(dat_period, env_vars_map, cutoff = 0.75)
  preds <- filtered$kept
  
  message("Running diagnostic period comparison: ", period_name)
  
  res <- run_daily_pls_sdm(
    data = dat_period,
    predictors = preds
  )
  
  all_results[[paste(period_name, "full", "pls", sep = "_")]] <- res
}

auc_summary <- map_dfr(names(all_results), function(nm) {
  
  res <- all_results[[nm]]
  if (is.null(res)) return(NULL)
  
  tibble(
    analysis = nm,
    auc = res$auc,
    youden_threshold = res$youden_threshold,
    n_predictors = length(res$predictors),
    predictors = paste(res$predictors, collapse = "; ")
  )
})

write.csv(
  auc_summary,
  file.path(out_dir, "Table_2_period_level_AUC_summary_diagnostic.csv"),
  row.names = FALSE
)


# 13) EFFORT DIAGNOSTIC — TOT
effort_vars <- c("n_detections_total", "n_species_detected", "mean_confidence")

base_preds_main <- filter_predictors(data_main, env_vars_map, cutoff = 0.75)$kept

effort_diag_main_pls <- run_daily_pls_sdm(
  data = data_main,
  predictors = c(base_preds_main, effort_vars)
)

effort_diag_table <- tibble(
  period = main_period,
  model = "pls",
  auc_with_effort_proxy = ifelse(
    is.null(effort_diag_main_pls),
    NA,
    effort_diag_main_pls$auc
  )
)

write.csv(
  effort_diag_table,
  file.path(out_dir, "effort_diagnostic_TOT_PLS_only.csv"),
  row.names = FALSE
)


# 14) AREA-SPECIFIC PLS MODELS — TOT
run_area_specific_pls <- function(data, area_name, period_name = "TOT",
                                  predictor_set = env_vars_map) {
  
  dat_area <- data %>%
    filter(biotope == area_name)
  
  if (length(unique(dat_area$id_area)) < 3) {
    message("Skipping ", area_name, ": fewer than 3 loggers.")
    return(NULL)
  }
  
  filtered <- filter_predictors(dat_area, predictor_set, cutoff = 0.75)
  preds <- filtered$kept
  
  message("Area model: ", area_name, " | ", period_name)
  
  res <- run_daily_pls_sdm(
    data = dat_area,
    predictors = preds
  )
  
  list(
    area = area_name,
    period = period_name,
    predictors = preds,
    removed_predictors = filtered$removed,
    model = res
  )
}

area_models_main_pls <- list()

for (area in names(area_rasters)) {
  
  area_models_main_pls[[area]] <- run_area_specific_pls(
    data = data_main,
    area_name = area,
    period_name = main_period,
    predictor_set = env_vars_map
  )
}


# 15) RASTER PREDICTION FUNCTION
predict_sdm_raster_pls <- function(raster_stack, model_object, predictors) {
  
  missing_layers <- setdiff(predictors, names(raster_stack))
  
  if (length(missing_layers) > 0) {
    stop("Missing raster layers: ", paste(missing_layers, collapse = ", "))
  }
  
  r <- raster_stack[[predictors]]
  names(r) <- predictors
  
  pred <- terra::predict(
    r,
    model_object,
    fun = function(model, data) {
      data <- as.data.frame(data)
      p <- predict(model, newdata = data, type = "prob")
      p[, "presence"]
    },
    na.rm = TRUE
  )
  
  names(pred) <- "suitability"
  return(pred)
}


# 16) SUITABILITY MAPS 
map_outputs <- list()

for (area in names(area_models_main_pls)) {
  
  area_obj <- area_models_main_pls[[area]]
  
  if (is.null(area_obj)) next
  if (is.null(area_obj$model)) next
  
  r <- area_rasters[[area]]
  
  pred_r <- predict_sdm_raster_pls(
    raster_stack = r,
    model_object = area_obj$model$model,
    predictors = area_obj$model$predictors
  )
  
  fname <- file.path(
    out_dir,
    paste0(main_period, "_full_", area, "_pls_suitability.tif")
  )
  
  writeRaster(pred_r, fname, overwrite = TRUE)
  
  map_outputs[[paste(main_period, "full", area, "pls", sep = "_")]] <- pred_r
}


# 17) EXTRAPOLATION / NOVELTY MAPS
make_extrapolation_map <- function(raster_stack, training_data, predictors) {
  
  r <- raster_stack[[predictors]]
  flags <- list()
  
  for (v in predictors) {
    
    v_min <- min(training_data[[v]], na.rm = TRUE)
    v_max <- max(training_data[[v]], na.rm = TRUE)
    
    flags[[v]] <- (r[[v]] < v_min) | (r[[v]] > v_max)
  }
  
  novelty <- Reduce(`+`, flags)
  names(novelty) <- "n_predictors_outside_training_range"
  
  return(novelty)
}

extrapolation_outputs <- list()

for (area in names(area_rasters)) {
  
  dat_area <- data_main %>% filter(biotope == area)
  
  if (nrow(dat_area) == 0) next
  if (is.null(area_models_main_pls[[area]])) next
  if (is.null(area_models_main_pls[[area]]$model)) next
  
  preds <- area_models_main_pls[[area]]$model$predictors
  
  novelty <- make_extrapolation_map(
    raster_stack = area_rasters[[area]],
    training_data = dat_area,
    predictors = preds
  )
  
  fname <- file.path(
    out_dir,
    paste0(main_period, "_full_", area, "_extrapolation_novelty.tif")
  )
  
  writeRaster(novelty, fname, overwrite = TRUE)
  extrapolation_outputs[[area]] <- novelty
}


# 18) AREA-SPECIFIC PERFORMANCE TABLE — TOT
final_auc_table <- map_dfr(names(area_models_main_pls), function(area) {
  
  obj <- area_models_main_pls[[area]]
  
  if (is.null(obj)) return(NULL)
  if (is.null(obj$model)) return(NULL)
  
  tibble(
    area = area,
    period = main_period,
    predictor_set = "full_resampled90_rao",
    model = "pls",
    auc = obj$model$auc,
    youden_threshold = obj$model$youden_threshold,
    n_loggers = length(unique(obj$model$data$id_area)),
    n_obs = nrow(obj$model$data),
    predictors = paste(obj$model$predictors, collapse = "; ")
  )
})

write.csv(
  final_auc_table,
  file.path(out_dir, "Table_3_area_specific_AUC_summary_TOT.csv"),
  row.names = FALSE
)


# 19) VARIABLE IMPORTANCE — TOT
importance_all <- map_dfr(names(area_models_main_pls), function(area) {
  
  obj <- area_models_main_pls[[area]]
  
  if (is.null(obj)) return(NULL)
  if (is.null(obj$model)) return(NULL)
  
  vi <- varImp(obj$model$model)$importance
  vi$variable <- rownames(vi)
  rownames(vi) <- NULL
  
  vi %>%
    mutate(
      area = area,
      predictor_set = "full_resampled90_rao",
      model = "pls",
      period = main_period
    )
})

write.csv(
  importance_all,
  file.path(out_dir, "variable_importance_area_specific_TOT_PLS_only.csv"),
  row.names = FALSE
)


# 20) DIAGNOSTIC SUMMARY
diagnostic_summary <- tibble(
  analysis = c(
    "site_level_leave_one_logger_out_PLS_TOT",
    "site_level_leave_one_area_out_PLS_TOT",
    "moran_residuals_site_level_PLS_TOT"
  ),
  statistic = c(
    "correlation",
    "correlation",
    "Moran_I"
  ),
  value = c(
    site_loocv_pls$cor,
    area_transfer_pls$cor,
    moran_site$observed
  ),
  rmse = c(
    site_loocv_pls$rmse,
    area_transfer_pls$rmse,
    NA
  ),
  p_value = c(
    NA,
    NA,
    moran_site$p_value
  )
)

write.csv(
  diagnostic_summary,
  file.path(out_dir, "Table_4_validation_diagnostics_TOT.csv"),
  row.names = FALSE
)


# 21) TABLES FOR PAPER
table_acoustic_summary <- logger_summary %>%
  mutate(
    prevalence = round(prevalence, 3),
    mean_detections_total = round(mean_detections_total, 1),
    mean_species_detected = round(mean_species_detected, 1)
  ) %>%
  arrange(biotope, id_area)

write.csv(
  table_acoustic_summary,
  file.path(out_dir, "Table_1_acoustic_sampling_summary.csv"),
  row.names = FALSE
)

table_period_auc <- auc_summary %>%
  mutate(
    auc = round(auc, 3),
    youden_threshold = round(youden_threshold, 3)
  )

write.csv(
  table_period_auc,
  file.path(out_dir, "Table_2_period_level_AUC_summary_diagnostic.csv"),
  row.names = FALSE
)

table_area_auc <- final_auc_table %>%
  mutate(
    auc = round(auc, 3),
    youden_threshold = round(youden_threshold, 3)
  ) %>%
  arrange(area)

write.csv(
  table_area_auc,
  file.path(out_dir, "Table_3_area_specific_AUC_summary_TOT.csv"),
  row.names = FALSE
)

table_diagnostics <- diagnostic_summary %>%
  mutate(
    value = round(value, 3),
    rmse = round(rmse, 3),
    p_value = signif(p_value, 3)
  )

write.csv(
  table_diagnostics,
  file.path(out_dir, "Table_4_validation_diagnostics_TOT.csv"),
  row.names = FALSE
)


# 22) FIGURES FOR PAPER
p_perf <- ggplot(final_auc_table, aes(x = area, y = auc, fill = area)) +
  geom_col(width = 0.7) +
  theme_bw() +
  labs(
    x = "Study area",
    y = "Leave-one-logger-out AUC",
    fill = "Area",
    title = "Area-specific PLS model performance using full-year data",
    subtitle = "LiDAR + Sentinel-2 + Rao predictors"
  ) + theme(
    axis.text.x = element_text(size = 12),  # valori asse X
    axis.text.y = element_text(size = 12),  # valori asse Y
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )

ggsave(
  file.path(out_dir, "Figure_3_area_specific_model_performance_TOT.png"),
  p_perf,
  width = 8,
  height = 5,
  dpi = 300
)

importance_plot_data <- importance_all %>%
  rename(importance = Overall) %>%
  group_by(area) %>%
  slice_max(order_by = importance, n = 10, with_ties = FALSE) %>%
  ungroup()

p_imp <- ggplot(
  importance_plot_data,
  aes(x = reorder(variable, importance), y = importance)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ area, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Predictor",
    y = "Variable importance",
    title = "Top predictors of full-year acoustic suitability",
    subtitle = "Area-specific PLS models, resampled 90 m predictors"
  ) + theme(
    strip.text = element_text(size = 14, face = "bold"), 
    axis.text.x = element_text(size = 10),  # valori asse X
    axis.text.y = element_text(size = 10)
  )

ggsave(
  file.path(out_dir, "Figure_4_variable_importance_area_specific_TOT.png"),
  p_imp,
  width = 11,
  height = 6,
  dpi = 300
)


# 23) SUITABILITY AND NOVELTY MAP PLOT
ahr_suit <- rast(file.path(out_dir, "TOT_full_ahr_pls_suitability.tif"))
slud_suit <- rast(file.path(out_dir, "TOT_full_slud_pls_suitability.tif"))
vals_suit <- rast(file.path(out_dir, "TOT_full_vals_pls_suitability.tif"))

ahr_ex <- rast(file.path(out_dir, "TOT_full_ahr_extrapolation_novelty.tif"))
slud_ex <- rast(file.path(out_dir, "TOT_full_slud_extrapolation_novelty.tif"))
vals_ex <- rast(file.path(out_dir, "TOT_full_vals_extrapolation_novelty.tif"))

png(
  file.path(out_dir, "Figure_5_suitability_maps_TOT.png"),
  width = 2200,
  height = 800,
  res = 200
)

par(mfrow = c(2, 2), mar = c(4, 4, 4, 5))
plot(ahr_suit,  main = "Ahr — suitability",  col = viridis(100), range = c(0, 1), cex.main = 1.5)
plot(slud_suit, main = "Slud — suitability", col = viridis(100), range = c(0, 1), cex.main = 1.5)
plot(vals_suit, main = "Vals — suitability", col = viridis(100), range = c(0, 1), cex.main = 1.5)

png(
  file.path(out_dir, "Figure_5_suitability_maps_TOT.png"),
  width = 2200,
  height = 800,
  res = 200
)


png(
  file.path(out_dir, "Figure_6_extrapolation_novelty_maps_TOT.png"),
  width = 2200,
  height = 800,
  res = 200
)

par(mfrow = c(1, 3), mar = c(4, 4, 4, 5))
plot(ahr_ex,  main = "Ahr — novelty")
plot(slud_ex, main = "Slud — novelty")
plot(vals_ex, main = "Vals — novelty")
dev.off()


# plot carini
par(
  mfrow = c(2, 2),
  mar = c(4, 4, 2, 5),   # margini singoli
  oma = c(0, 0, 4, 0)    # margine esterno sopra
)

plot(
  ahr_suit,
  main = "Ahr",
  col = viridis(100),
  range = c(0, 1),
  cex.main = 1
)

plot(
  slud_suit,
  main = "Slud",
  col = viridis(100),
  range = c(0, 1),
  cex.main = 1
)

plot(
  vals_suit,
  main = "Vals",
  col = viridis(100),
  range = c(0, 1),
  cex.main = 1
)

mtext(
  "Predicted habitat suitability for Picus viridis",
  outer = TRUE,
  cex = 1.2,
  font = 2,
  line = 1
)

# novelty
library(viridis)

max_novelty <- max(
  values(ahr_ex),
  values(slud_ex),
  values(vals_ex),
  na.rm = TRUE
)

cols_novelty <- viridis(max_novelty + 1)

par(
  mfrow = c(1, 3),
  mar = c(4, 4, 2, 2),
  oma = c(0, 0, 4, 4)
)

plot(
  ahr_ex,
  main = "Ahr",
  col = cols_novelty,
  breaks = seq(-0.5, max_novelty + 0.5, by = 1),
  legend = FALSE,
  cex.main = 1
)

plot(
  slud_ex,
  main = "Slud",
  col = cols_novelty,
  breaks = seq(-0.5, max_novelty + 0.5, by = 1),
  legend = FALSE,
  cex.main = 1
)

plot(
  vals_ex,
  main = "Vals",
  col = cols_novelty,
  breaks = seq(-0.5, max_novelty + 0.5, by = 1),
  cex.main = 1, 
  plg = list(
    cex = 1.3
  )
)

mtext(
  "Environmental novelty: number of predictors outside the training range",
  outer = TRUE,
  cex = 1.2,
  font = 2,
  line = 1
)


# 24) FINAL CHECK
print("=== TABLE 1: Acoustic summary ===")
print(table_acoustic_summary, n = Inf)

print("=== TABLE 2: Period AUC diagnostic ===")
print(table_period_auc)

print("=== TABLE 3: Area-specific AUC TOT ===")
print(table_area_auc)

print("=== TABLE 4: Diagnostics TOT ===")
print(table_diagnostics)

print("=== Moran site-level residuals ===")
print(moran_site)

message("Pipeline completed.")
message("Outputs saved in: ", out_dir)

list.files(out_dir, pattern = ".tif")

ahr_suit <- rast("outputs_revision_PLS_TOT_resampled90_rao/TOT_full_ahr_pls_suitability.tif")
slud_suit <- rast("outputs_revision_PLS_TOT_resampled90_rao/TOT_full_slud_pls_suitability.tif")
vals_suit <- rast("outputs_revision_PLS_TOT_resampled90_rao/TOT_full_vals_pls_suitability.tif")

ahr_ex <- rast("outputs_revision_PLS_TOT_resampled90_rao/TOT_full_ahr_extrapolation_novelty.tif")
slud_ex <- rast("outputs_revision_PLS_TOT_resampled90_rao/TOT_full_slud_extrapolation_novelty.tif")
vals_ex <- rast("outputs_revision_PLS_TOT_resampled90_rao/TOT_full_vals_extrapolation_novelty.tif")

par(mfrow = c(1,3))
plot(ahr_suit, main = "ahr suitability")
plot(slud_suit, main = "slud suitability")
plot(vals_suit, main = "vals suitability")

par(mfrow = c(1,3))
plot(ahr_ex, main = "ahr extrap. novelty")
plot(slud_ex, main = "slud extrap. novelty")
plot(vals_ex, main = "vals extrap. novelty")


# test
library(terra)
library(ggplot2)
library(tidyterra)
library(ggspatial)
library(patchwork)
library(scales)

make_map <- function(r, title, scale_width = 0.25) {
  ggplot() +
    geom_spatraster(data = r, na.rm = TRUE) +
    scale_fill_viridis_c(
      name = "Suitability",
      labels = label_number(accuracy = 0.01),
      na.value = "transparent"
    ) +
    annotation_scale(
      location = "bl",
      width_hint = scale_width,
      text_cex = 0.7
    ) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      height = unit(0.6, "cm"),
      width  = unit(0.6, "cm"),
      pad_x = unit(0.15, "cm"),
      pad_y = unit(0.15, "cm"),
      style = north_arrow_minimal
    ) +
    labs(title = title) +
    coord_sf(
      crs = terra::crs(r),
      datum = terra::crs(r)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.title = element_blank(),
      axis.text.x = element_text(
        angle = 45,      # Rotation angle in degrees
        hjust = 1,       # Horizontal justification
        vjust = 1        # Vertical justification
      )
    )
}

p1 <- make_map(ahr_suit,  "Greti dell'Aurino suitability",  0.35)
p2 <- make_map(slud_suit, "Sludern suitability", 0.25)
p3 <- make_map(vals_suit, "Valsura Delta suitability", 0.25)

final_plot <- (p1 + p2) / (p3 + plot_spacer())

ggsave(
  "suitability_maps.png",
  plot = final_plot,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "transparent"
)

# --- 
library(xtable)
library(knitr)
library(gt)

print(
  xtable(table_acoustic_summary, caption = "Acoustic summary", digits = 2),
  include.rownames = FALSE
)

kable(table_period_auc, caption = "Period AUC diagnosticy", digits = 2)
kable(table_area_auc, caption = "Area-specific AUC", digits = 2)

table_acoustic_summary %>%
  gt() %>%
  tab_header(
    title = "Acoustic summary"
  )


# riscrivi tabelle

list.files(out_dir, pattern = 'csv')

write.csv2(
  table_acoustic_summary,
  file.path(out_dir, "Table_1_acoustic_sampling_summary_excel.csv"),
  row.names = FALSE
)

write.csv2(
  table_period_auc,
  file.path(out_dir, "Table_2_period_level_AUC_summary_diagnostic_exc.csv"),
  row.names = FALSE
)

write.csv2(
  table_area_auc,
  file.path(out_dir, "Table_3_area_specific_AUC_summary_TOT_exc.csv"),
  row.names = FALSE
)

write.csv2(
  table_diagnostics,
  file.path(out_dir, "Table_4_validation_diagnostics_TOT.csv"),
  row.names = FALSE
)

# moran_site
write.csv2(
  moran_site,
  file.path(out_dir, "moran_residuals_site_level_LOOCV_PLS_TOT_exc.csv"),
  row.names = FALSE
)

library(dplyr)
library(tidyr)

# dataset con una riga per logger + coordinate + predictors
logger_env <- site_data %>%
  select(id_area, biotope, x, y, all_of(env_vars)) %>%
  distinct()

# standardizzo predictor ambientali
env_scaled <- logger_env %>%
  select(all_of(env_vars)) %>%
  scale()

# distanze ambientali
env_dist <- as.matrix(dist(env_scaled))

# distanze spaziali
coords <- logger_env %>%
  select(x, y) %>%
  as.matrix()

spatial_dist <- as.matrix(dist(coords))

# logger esclusi
removed_loggers <- c("vals_moth9", "vals_moth4")

redundancy_table <- lapply(removed_loggers, function(id) {
  
  i <- which(logger_env$id_area == id)
  
  candidate_idx <- which(
    logger_env$biotope == logger_env$biotope[i] &
      logger_env$id_area != id &
      !logger_env$id_area %in% removed_loggers
  )
  
  nearest_idx <- candidate_idx[which.min(spatial_dist[i, candidate_idx])]
  most_similar_idx <- candidate_idx[which.min(env_dist[i, candidate_idx])]
  
  data.frame(
    removed_logger = id,
    biotope = logger_env$biotope[i],
    nearest_retained_logger = logger_env$id_area[nearest_idx],
    distance_to_nearest_m = round(spatial_dist[i, nearest_idx], 1),
    most_environmentally_similar_logger = logger_env$id_area[most_similar_idx],
    environmental_distance = round(env_dist[i, most_similar_idx], 3),
    reason = "Spatially close and environmentally redundant"
  )
}) %>%
  bind_rows()

redundancy_table

write.csv2(
  redundancy_table,
  file.path(out_dir, "Table_S1_removed_logger_redundancy.csv"),
  row.names = FALSE
)
