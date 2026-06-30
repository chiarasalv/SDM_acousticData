options(warn = -1)
# ============================================================
# Reproducible PLSR workflow — Picus viridis
# Daily processed dataset + raster predictors
# Import data from Zenodo!
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

# 1. Paths and input files
# SDM_acousticData/
# │
# ├── README.md
# ├── LICENSE
# ├── scripts/
# ├── outputs/
# ├── zenodo_data/
# │      PicusViridis_daily_model_data.csv
# │      logger_coordinates.csv
# │      ahr_metrics_resampled90_masked.tif
# │      slud_metrics_resampled90_masked.tif
# │      vals_metrics_resampled90_masked.tif
# │      README_DATA.txt

project_dir  <- getwd()
prepared_dir <- file.path(project_dir, "zenodo_data")
out_dir      <- file.path(project_dir, "output")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_files <- c(
  "PicusViridis_daily_model_data.csv",
  "logger_coordinates.csv",
  "ahr_metrics_resampled90_masked.tif",
  "slud_metrics_resampled90_masked.tif",
  "vals_metrics_resampled90_masked.tif"
)

missing_files <- required_files[
  !file.exists(file.path(prepared_dir, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "Missing required input files in: ", prepared_dir, "\n",
    paste(missing_files, collapse = "\n")
  )
}

# Load inputs
model_data_raw <- read.csv(
  file.path(prepared_dir, "PicusViridis_daily_model_data.csv"),
  stringsAsFactors = FALSE
)

logger_points_raw <- read.csv(
  file.path(prepared_dir, "logger_coordinates.csv"),
  stringsAsFactors = FALSE
)

area_rasters <- list(
  ahr  = rast(file.path(prepared_dir, "ahr_metrics_resampled90_masked.tif")),
  slud = rast(file.path(prepared_dir, "slud_metrics_resampled90_masked.tif")),
  vals = rast(file.path(prepared_dir, "vals_metrics_resampled90_masked.tif"))
)

# 2. Settings
target_species <- "Picus viridis"
main_period <- "TOT"

loggers_to_remove <- c("vals_moth4", "vals_moth9")

env_vars <- c(
  "vegetation_cover", "canopy_volume", "chm_rugosity", "num_trees",
  "zmax", "zmean", "zsd", "zskew", "zkurt",
  "pzabovezmean", "pzabove2",
  paste0("zq", seq(5, 95, 5)),
  paste0("zpcum", 1:9),
  "NDVI_mean", "NDVI_sd", "NDVI_max",
  "rao_chm", "rao_ndvi"
)


# 3. Input checks and formatting
required_columns <- c(
  "id_area", "biotope", "logger", "Date", "period",
  "target_species", "presence",
  "n_detections_total", "n_species_detected", "mean_confidence",
  env_vars
)

missing_columns <- setdiff(required_columns, names(model_data_raw))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns in PicusViridis_daily_model_data.csv:\n",
    paste(missing_columns, collapse = "\n")
  )
}

model_data_all <- model_data_raw %>%
  mutate(
    Date = as.Date(Date),
    id_area = as.character(id_area),
    biotope = as.character(biotope),
    logger = as.character(logger),
    period = as.character(period),
    target_species = as.character(target_species),
    presence = as.integer(presence)
  )

if (!all(model_data_all$target_species == target_species)) {
  stop("The input dataset contains target_species values different from Picus viridis.")
}

if (!all(model_data_all$presence %in% c(0, 1))) {
  stop("The presence column must contain only 0/1 values.")
}

if (!all(c("id_area", "x", "y") %in% names(logger_points_raw))) {
  stop("logger_coordinates.csv must contain at least: id_area, x, y.")
}


# 4. Preliminary redundancy table BEFORE removal
build_site_dataset <- function(model_data, logger_points, env_vars) {
  model_data %>%
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
}

site_data_pre_removal <- build_site_dataset(
  model_data = model_data_all,
  logger_points = logger_points_raw,
  env_vars = env_vars
)

vals_logger <- site_data_pre_removal %>%
  filter(biotope == "vals")

if (nrow(vals_logger) > 1) {
  
  coords <- vals_logger %>%
    select(x, y) %>%
    as.matrix()
  
  spatial_dist <- as.matrix(dist(coords))
  rownames(spatial_dist) <- vals_logger$id_area
  colnames(spatial_dist) <- vals_logger$id_area
  diag(spatial_dist) <- NA
  
  env_scaled <- vals_logger %>%
    select(all_of(env_vars)) %>%
    scale()
  
  env_dist <- as.matrix(dist(env_scaled))
  rownames(env_dist) <- vals_logger$id_area
  colnames(env_dist) <- vals_logger$id_area
  diag(env_dist) <- NA
  
  vals_redundancy_summary <- map_dfr(vals_logger$id_area, function(id) {
    
    nearest_spatial <- names(which.min(spatial_dist[id, ]))
    nearest_env <- names(which.min(env_dist[id, ]))
    
    tibble(
      id_area = id,
      decision = if_else(id %in% loggers_to_remove, "removed", "retained"),
      prevalence = vals_logger$presence_rate[vals_logger$id_area == id],
      mean_detections_total = vals_logger$mean_detections_total[vals_logger$id_area == id],
      mean_species_detected = vals_logger$mean_species_detected[vals_logger$id_area == id],
      nearest_spatial_logger = nearest_spatial,
      spatial_distance_m = round(spatial_dist[id, nearest_spatial], 1),
      nearest_environmental_logger = nearest_env,
      environmental_distance = round(env_dist[id, nearest_env], 3),
      redundancy_group = case_when(
        id %in% c("vals_moth3", "vals_moth4") ~ "G1",
        id %in% c("vals_moth8", "vals_moth9") ~ "G2",
        TRUE ~ NA_character_
      ),
      reason = if_else(
        id %in% loggers_to_remove,
        "Spatially and environmentally redundant within Valsura Delta",
        "Retained for spatial and environmental representation"
      )
    ) %>%
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
      )
  })
  
  write.csv(
    vals_redundancy_summary,
    file.path(out_dir, "Table_S1_valsura_logger_redundancy_assessment.csv"),
    row.names = FALSE
  )
}


# 5. Apply logger removal BEFORE all final datasets
model_data <- model_data_all %>%
  filter(!id_area %in% loggers_to_remove)

logger_points <- logger_points_raw %>%
  filter(!id_area %in% loggers_to_remove)

data_TOT <- model_data
data_PR  <- model_data %>% filter(period == "PR")
data_NR  <- model_data %>% filter(period == "NR")
data_main <- data_TOT

site_data <- build_site_dataset(
  model_data = data_main,
  logger_points = logger_points,
  env_vars = env_vars
)

logger_summary <- data_main %>%
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

check_removed_absent <- function(data, name) {
  found <- intersect(unique(data$id_area), loggers_to_remove)
  if (length(found) > 0) {
    stop(
      "Removed loggers found in ", name, ": ",
      paste(found, collapse = ", ")
    )
  }
}

check_removed_absent(model_data, "model_data")
check_removed_absent(data_TOT, "data_TOT")
check_removed_absent(data_PR, "data_PR")
check_removed_absent(data_NR, "data_NR")
check_removed_absent(site_data, "site_data")
check_removed_absent(logger_summary, "logger_summary")

message("Check passed: removed loggers are absent from all final datasets.")


# 6. Helper functions
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

run_daily_plsr_sdm <- function(data, predictors) {
  
  dat <- data %>%
    select(id_area, biotope, presence, all_of(predictors)) %>%
    drop_na()
  
  dat$presence <- factor(
    dat$presence,
    levels = c(0, 1),
    labels = c("absence", "presence")
  )
  
  if (length(unique(dat$presence)) < 2) {
    message("Skipping PLSR: only one response class.")
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
    savePredictions = "final"
  )
  
  tune_grid <- expand.grid(
    ncomp = 1:min(5, length(predictors))
  )
  
  set.seed(123)
  
  mod <- train(
    presence ~ .,
    data = dat %>% select(-id_area, -biotope),
    method = "pls",
    metric = "ROC",
    trControl = ctrl,
    tuneGrid = tune_grid
  )
  
  pred <- mod$pred %>%
    filter(!is.na(presence))
  
  if (length(unique(pred$obs)) < 2) {
    message("Skipping PLSR: only one class in predictions.")
    return(NULL)
  }
  
  roc_obj <- pROC::roc(
    response = pred$obs,
    predictor = pred$presence,
    levels = c("absence", "presence"),
    quiet = TRUE
  )
  
  best_thr <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = "threshold"
  )
  
  list(
    model = mod,
    predictions = pred,
    auc = as.numeric(pROC::auc(roc_obj)),
    youden_threshold = as.numeric(best_thr),
    predictors = predictors,
    data = dat,
    best_ncomp = mod$bestTune$ncomp,
    tuning_results = mod$results
  )
}

run_site_loocv_plsr <- function(site_data, predictors) {
  
  dat <- site_data %>%
    select(id_area, biotope, presence_rate, x, y, all_of(predictors)) %>%
    drop_na()
  
  preds <- rep(NA_real_, nrow(dat))
  
  for (i in seq_len(nrow(dat))) {
    
    train_dat <- dat[-i, ]
    test_dat  <- dat[i, ]
    
    ncomp_i <- min(3, length(predictors), nrow(train_dat) - 1)
    if (ncomp_i < 1) next
    
    set.seed(123)
    
    mod <- train(
      presence_rate ~ .,
      data = train_dat %>% select(-id_area, -biotope, -x, -y),
      method = "pls",
      trControl = trainControl(method = "none"),
      tuneGrid = data.frame(ncomp = ncomp_i)
    )
    
    preds[i] <- predict(mod, newdata = test_dat)
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

run_leave_one_area_plsr <- function(site_data, predictors) {
  
  dat <- site_data %>%
    select(id_area, biotope, presence_rate, all_of(predictors)) %>%
    drop_na()
  
  all_preds <- map_dfr(unique(dat$biotope), function(area) {
    
    train_dat <- dat %>% filter(biotope != area)
    test_dat  <- dat %>% filter(biotope == area)
    
    ncomp_i <- min(3, length(predictors), nrow(train_dat) - 1)
    if (ncomp_i < 1) return(NULL)
    
    set.seed(123)
    
    mod <- train(
      presence_rate ~ .,
      data = train_dat %>% select(-id_area, -biotope),
      method = "pls",
      trControl = trainControl(method = "none"),
      tuneGrid = data.frame(ncomp = ncomp_i)
    )
    
    test_dat %>%
      mutate(
        predicted = predict(mod, newdata = test_dat),
        residual = presence_rate - predicted,
        left_out_area = area
      )
  })
  
  list(
    predictions = all_preds,
    cor = cor(all_preds$presence_rate, all_preds$predicted, use = "complete.obs"),
    rmse = sqrt(mean((all_preds$presence_rate - all_preds$predicted)^2, na.rm = TRUE))
  )
}

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

predict_sdm_raster_plsr <- function(raster_stack, model_object, predictors) {
  
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
      p <- predict(model, newdata = as.data.frame(data), type = "prob")
      p[, "presence"]
    },
    na.rm = TRUE
  )
  
  names(pred) <- "relative_acoustic_habitat_suitability"
  pred
}

make_novelty_map <- function(raster_stack, training_data, predictors) {
  
  r <- raster_stack[[predictors]]
  
  flags <- map(predictors, function(v) {
    v_min <- min(training_data[[v]], na.rm = TRUE)
    v_max <- max(training_data[[v]], na.rm = TRUE)
    (r[[v]] < v_min) | (r[[v]] > v_max)
  })
  
  novelty <- Reduce(`+`, flags)
  names(novelty) <- "n_predictors_outside_training_range"
  
  novelty
}


# 7. Tables: acoustic summary
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


# 8. Period comparison: diagnostic only
period_datasets <- list(
  TOT = data_TOT,
  PR  = data_PR,
  NR  = data_NR
)

period_results <- map(period_datasets, function(dat_period) {
  
  filtered <- filter_predictors(dat_period, env_vars, cutoff = 0.75)
  
  run_daily_plsr_sdm(
    data = dat_period,
    predictors = filtered$kept
  )
})

period_auc_table <- map_dfr(names(period_results), function(period_name) {
  
  res <- period_results[[period_name]]
  if (is.null(res)) return(NULL)
  
  tibble(
    analysis = paste(period_name, "full", "PLSR", sep = "_"),
    period = period_name,
    auc = round(res$auc, 3),
    youden_threshold = round(res$youden_threshold, 3),
    n_predictors = length(res$predictors),
    predictors = paste(res$predictors, collapse = "; ")
  )
})

write.csv(
  period_auc_table,
  file.path(out_dir, "Table_2_period_level_AUC_summary_diagnostic.csv"),
  row.names = FALSE
)


# 9. Site-level validation and transferability — TOT
site_predictors <- filter_predictors(site_data, env_vars, cutoff = 0.75)$kept

site_loocv_plsr <- run_site_loocv_plsr(
  site_data = site_data,
  predictors = site_predictors
)

area_transfer_plsr <- run_leave_one_area_plsr(
  site_data = site_data,
  predictors = site_predictors
)

moran_site <- run_moran_residuals(site_loocv_plsr$predictions)

write.csv(
  site_loocv_plsr$predictions,
  file.path(out_dir, "site_level_LOLO_predictions_PLSR_TOT.csv"),
  row.names = FALSE
)

write.csv(
  area_transfer_plsr$predictions,
  file.path(out_dir, "leave_one_area_out_predictions_PLSR_TOT.csv"),
  row.names = FALSE
)

diagnostic_summary <- tibble(
  analysis = c(
    "site_level_leave_one_logger_out_PLSR_TOT",
    "site_level_leave_one_area_out_PLSR_TOT",
    "moran_residuals_site_level_PLSR_TOT"
  ),
  statistic = c("correlation", "correlation", "Moran_I"),
  value = c(site_loocv_plsr$cor, area_transfer_plsr$cor, moran_site$observed),
  rmse = c(site_loocv_plsr$rmse, area_transfer_plsr$rmse, NA),
  p_value = c(NA, NA, moran_site$p_value)
) %>%
  mutate(
    value = round(value, 3),
    rmse = round(rmse, 3),
    p_value = signif(p_value, 3)
  )

write.csv(
  diagnostic_summary,
  file.path(out_dir, "Table_4_validation_diagnostics_TOT.csv"),
  row.names = FALSE
)


# 10. Area-specific PLSR models
run_area_specific_plsr <- function(data, area_name, predictor_set = env_vars) {
  
  dat_area <- data %>%
    filter(biotope == area_name)
  
  if (length(unique(dat_area$id_area)) < 3) {
    message("Skipping ", area_name, ": fewer than 3 loggers.")
    return(NULL)
  }
  
  filtered <- filter_predictors(dat_area, predictor_set, cutoff = 0.75)
  
  res <- run_daily_plsr_sdm(
    data = dat_area,
    predictors = filtered$kept
  )
  
  list(
    area = area_name,
    predictors = filtered$kept,
    removed_predictors = filtered$removed,
    model = res
  )
}

area_models <- map(
  names(area_rasters),
  ~ run_area_specific_plsr(data_main, .x, env_vars)
)

names(area_models) <- names(area_rasters)

area_auc_table <- map_dfr(names(area_models), function(area) {
  
  obj <- area_models[[area]]
  if (is.null(obj) || is.null(obj$model)) return(NULL)
  
  tibble(
    area = area,
    period = main_period,
    model = "PLSR",
    auc = round(obj$model$auc, 3),
    youden_threshold = round(obj$model$youden_threshold, 3),
    best_ncomp = obj$model$best_ncomp,
    n_loggers = length(unique(obj$model$data$id_area)),
    n_obs = nrow(obj$model$data),
    predictors = paste(obj$model$predictors, collapse = "; ")
  )
})

write.csv(
  area_auc_table,
  file.path(out_dir, "Table_3_area_specific_AUC_summary_TOT.csv"),
  row.names = FALSE
)

component_tuning_table <- map_dfr(names(area_models), function(area) {
  
  obj <- area_models[[area]]
  if (is.null(obj) || is.null(obj$model)) return(NULL)
  
  obj$model$tuning_results %>%
    mutate(
      area = area,
      best_ncomp = obj$model$best_ncomp
    )
})

write.csv(
  component_tuning_table,
  file.path(out_dir, "Table_S3_PLSR_component_tuning.csv"),
  row.names = FALSE
)

importance_table <- map_dfr(names(area_models), function(area) {
  
  obj <- area_models[[area]]
  if (is.null(obj) || is.null(obj$model)) return(NULL)
  
  vi <- varImp(obj$model$model)$importance
  vi$predictor <- rownames(vi)
  rownames(vi) <- NULL
  
  vi %>%
    mutate(
      area = area,
      period = main_period,
      model = "PLSR"
    )
})

write.csv(
  importance_table,
  file.path(out_dir, "variable_importance_area_specific_TOT_PLSR.csv"),
  row.names = FALSE
)


# 11. Raster prediction, novelty maps, and range diagnostics
suitability_outputs <- list()
novelty_outputs <- list()
range_diagnostics <- list()

for (area in names(area_models)) {
  
  obj <- area_models[[area]]
  if (is.null(obj) || is.null(obj$model)) next
  
  dat_area <- data_main %>% filter(biotope == area)
  preds <- obj$model$predictors
  
  suitability <- predict_sdm_raster_plsr(
    raster_stack = area_rasters[[area]],
    model_object = obj$model$model,
    predictors = preds
  )
  
  novelty <- make_novelty_map(
    raster_stack = area_rasters[[area]],
    training_data = dat_area,
    predictors = preds
  )
  
  writeRaster(
    suitability,
    file.path(out_dir, paste0("TOT_", area, "_PLSR_suitability.tif")),
    overwrite = TRUE
  )
  
  writeRaster(
    novelty,
    file.path(out_dir, paste0("TOT_", area, "_PLSR_environmental_novelty.tif")),
    overwrite = TRUE
  )
  
  suitability_outputs[[area]] <- suitability
  novelty_outputs[[area]] <- novelty
  
  range_diagnostics[[area]] <- map_dfr(preds, function(v) {
    
    raster_values <- values(area_rasters[[area]][[v]], na.rm = TRUE)
    train_values <- dat_area[[v]]
    
    training_min <- min(train_values, na.rm = TRUE)
    training_max <- max(train_values, na.rm = TRUE)
    
    n_below <- sum(raster_values < training_min, na.rm = TRUE)
    n_above <- sum(raster_values > training_max, na.rm = TRUE)
    n_cells <- length(raster_values)
    
    tibble(
      area = area,
      predictor = v,
      training_min = training_min,
      training_max = training_max,
      raster_min = min(raster_values, na.rm = TRUE),
      raster_max = max(raster_values, na.rm = TRUE),
      n_cells = n_cells,
      n_below = n_below,
      n_above = n_above,
      n_outside = n_below + n_above,
      pct_inside = round(100 * (1 - (n_below + n_above) / n_cells), 1),
      pct_outside = round(100 * (n_below + n_above) / n_cells, 1)
    )
  })
}

range_checks_df <- bind_rows(range_diagnostics)

write.csv(
  range_checks_df,
  file.path(out_dir, "Table_S2_predictor_range_diagnostics.csv"),
  row.names = FALSE
)

# 12. Figures
p_site_loocv <- ggplot(
  site_loocv_plsr$predictions,
  aes(x = presence_rate, y = predicted, color = biotope)
) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    x = "Observed full-year acoustic occurrence rate",
    y = "Predicted acoustic occurrence rate",
    color = "Study area"
  )

ggsave(
  file.path(out_dir, "Figure_3_site_level_LOLO_PLSR_TOT.png"),
  p_site_loocv,
  width = 7,
  height = 5,
  dpi = 300
)

p_area_transfer <- ggplot(
  area_transfer_plsr$predictions,
  aes(x = presence_rate, y = predicted, color = left_out_area)
) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    x = "Observed full-year acoustic occurrence rate",
    y = "Predicted acoustic occurrence rate",
    color = "Left-out area"
  )

ggsave(
  file.path(out_dir, "Figure_4_leave_one_area_out_PLSR_TOT.png"),
  p_area_transfer,
  width = 7,
  height = 5,
  dpi = 300
)

p_perf <- ggplot(area_auc_table, aes(x = area, y = auc, fill = area)) +
  geom_col(width = 0.7) +
  theme_bw() +
  labs(
    x = "Study area",
    y = "Leave-one-logger-out AUC",
    fill = "Study area"
  )

ggsave(
  file.path(out_dir, "Figure_5_area_specific_model_performance_TOT.png"),
  p_perf,
  width = 8,
  height = 5,
  dpi = 300
)

importance_plot_data <- importance_table %>%
  rename(importance = Overall) %>%
  group_by(area) %>%
  slice_max(order_by = importance, n = 10, with_ties = FALSE) %>%
  ungroup()

p_imp <- ggplot(
  importance_plot_data,
  aes(x = reorder(predictor, importance), y = importance)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ area, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Predictor",
    y = "Relative variable importance"
  )

ggsave(
  file.path(out_dir, "Figure_6_variable_importance_area_specific_TOT.png"),
  p_imp,
  width = 11,
  height = 6,
  dpi = 300
)

png(
  file.path(out_dir, "Figure_7_suitability_maps_TOT.png"),
  width = 2200,
  height = 800,
  res = 200
)

par(mfrow = c(1, 3), mar = c(4, 4, 3, 5))

plot(
  suitability_outputs$ahr,
  main = "Greti dell'Aurino (ahr)",
  col = viridis(100)
)

plot(
  suitability_outputs$slud,
  main = "Sluderno Alderwood (slud)",
  col = viridis(100)
)

plot(
  suitability_outputs$vals,
  main = "Valsura Delta (vals)",
  col = viridis(100)
)

dev.off()

max_novelty <- max(
  values(novelty_outputs$ahr),
  values(novelty_outputs$slud),
  values(novelty_outputs$vals),
  na.rm = TRUE
)

png(
  file.path(out_dir, "Figure_8_environmental_novelty_maps_TOT.png"),
  width = 2200,
  height = 800,
  res = 200
)

par(mfrow = c(1, 3), mar = c(4, 4, 3, 5))

plot(
  novelty_outputs$ahr,
  main = "Greti dell'Aurino (ahr)",
  col = viridis(max_novelty + 1),
  breaks = seq(-0.5, max_novelty + 0.5, by = 1)
)

plot(
  novelty_outputs$slud,
  main = "Sluderno Alderwood (slud)",
  col = viridis(max_novelty + 1),
  breaks = seq(-0.5, max_novelty + 0.5, by = 1)
)

plot(
  novelty_outputs$vals,
  main = "Valsura Delta (vals)",
  col = viridis(max_novelty + 1),
  breaks = seq(-0.5, max_novelty + 0.5, by = 1)
)

dev.off()


# 13. Final reproducibility checks
final_outputs <- c(
  "Table_1_acoustic_sampling_summary.csv",
  "Table_2_period_level_AUC_summary_diagnostic.csv",
  "Table_3_area_specific_AUC_summary_TOT.csv",
  "Table_4_validation_diagnostics_TOT.csv",
  "Table_S1_valsura_logger_redundancy_assessment.csv",
  "Table_S2_predictor_range_diagnostics.csv",
  "Table_S3_PLSR_component_tuning.csv",
  "Figure_3_site_level_LOLO_PLSR_TOT.png",
  "Figure_4_leave_one_area_out_PLSR_TOT.png",
  "Figure_5_area_specific_model_performance_TOT.png",
  "Figure_6_variable_importance_area_specific_TOT.png",
  "Figure_7_suitability_maps_TOT.png",
  "Figure_8_environmental_novelty_maps_TOT.png"
)

missing_outputs <- final_outputs[
  !file.exists(file.path(out_dir, final_outputs))
]

if (length(missing_outputs) > 0) {
  stop(
    "Some expected outputs were not generated:\n",
    paste(missing_outputs, collapse = "\n")
  )
}

message("Pipeline completed successfully.")
message("Outputs saved in: ", out_dir)
