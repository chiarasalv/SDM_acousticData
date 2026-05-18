library(dplyr)
library(tidyr)
library(terra)
library(caret)
library(pROC)
library(pls)
library(ggplot2)
library(purrr)
library(ape)
library(tibble)
library(viridis)

set.seed(123)

config_dir <- "config"
data_dir <- "data"
output_dir <- "outputs/pls_sdm"
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
raster_dir <- file.path(output_dir, "rasters")
environment_dir <- "environment"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(environment_dir, recursive = TRUE, showWarnings = FALSE)

analysis_config <- read.csv(file.path(config_dir, "analysis_config.csv"))
predictor_config <- read.csv(file.path(config_dir, "predictor_variables.csv"))
raster_inventory <- read.csv(file.path(config_dir, "raster_inventory.csv"))
site_exclusion <- read.csv(file.path(config_dir, "site_exclusion.csv"))

required_config_columns <- c("parameter", "value")
if (!all(required_config_columns %in% names(analysis_config))) {
  stop("analysis_config.csv must contain columns: parameter, value")
}

get_config_value <- function(config, parameter_name) {
  value <- config$value[config$parameter == parameter_name]
  if (length(value) == 0) stop("Missing parameter in analysis_config.csv: ", parameter_name)
  value[1]
}

input_table <- get_config_value(analysis_config, "input_table")
site_table <- get_config_value(analysis_config, "site_table")
target_species <- get_config_value(analysis_config, "target_species")
main_period <- get_config_value(analysis_config, "main_period")
breeding_start <- as.Date(get_config_value(analysis_config, "breeding_start"))
breeding_end <- as.Date(get_config_value(analysis_config, "breeding_end"))
correlation_cutoff <- as.numeric(get_config_value(analysis_config, "correlation_cutoff"))
minimum_loggers_per_area <- as.numeric(get_config_value(analysis_config, "minimum_loggers_per_area"))

required_predictor_columns <- "predictor"
if (!all(required_predictor_columns %in% names(predictor_config))) {
  stop("predictor_variables.csv must contain a column named predictor")
}

predictor_vars <- predictor_config$predictor

detection_data <- read.csv(file.path(data_dir, input_table))
site_data_raw <- read.csv(file.path(data_dir, site_table))

required_detection_columns <- c("site_id", "area_id", "date", "species", "confidence")
missing_detection_columns <- setdiff(required_detection_columns, names(detection_data))
if (length(missing_detection_columns) > 0) {
  stop("Input table is missing columns: ", paste(missing_detection_columns, collapse = ", "))
}

required_site_columns <- c("site_id", "area_id", "x", "y")
missing_site_columns <- setdiff(required_site_columns, names(site_data_raw))
if (length(missing_site_columns) > 0) {
  stop("Site table is missing columns: ", paste(missing_site_columns, collapse = ", "))
}

missing_predictors <- setdiff(predictor_vars, names(detection_data))
if (length(missing_predictors) > 0) {
  stop("Input table is missing predictor columns: ", paste(missing_predictors, collapse = ", "))
}

required_raster_columns <- c("area_id", "raster_path")
missing_raster_columns <- setdiff(required_raster_columns, names(raster_inventory))
if (length(missing_raster_columns) > 0) {
  stop("raster_inventory.csv is missing columns: ", paste(missing_raster_columns, collapse = ", "))
}

area_rasters <- raster_inventory %>%
  mutate(raster = map(file.path(data_dir, raster_path), rast)) %>%
  pull(raster, name = area_id)

writeLines(capture.output(sessionInfo()), file.path(environment_dir, "session_info.txt"))

prepare_detection_data <- function(data, breeding_start, breeding_end) {
  data %>%
    mutate(
      date = as.Date(date),
      period = case_when(
        date >= breeding_start & date <= breeding_end ~ "PR",
        TRUE ~ "NR"
      )
    )
}

build_daily_dataset <- function(data, target_species, predictor_vars) {
  daily_data <- data %>%
    group_by(site_id, area_id, date, period) %>%
    summarise(
      presence = as.integer(any(species == target_species)),
      n_detections_total = n(),
      n_species_detected = n_distinct(species),
      mean_confidence = mean(confidence, na.rm = TRUE),
      .groups = "drop"
    )
  
  site_predictors <- data %>%
    select(site_id, all_of(predictor_vars)) %>%
    distinct()
  
  predictor_check <- site_predictors %>%
    count(site_id) %>%
    filter(n > 1)
  
  if (nrow(predictor_check) > 0) {
    warning("Some sites have more than one row of predictor values. Keeping distinct combinations may duplicate site-day records.")
  }
  
  daily_data %>%
    left_join(site_predictors, by = "site_id") %>%
    drop_na()
}

build_site_dataset <- function(daily_data, site_data, predictor_vars) {
  daily_data %>%
    group_by(site_id, area_id) %>%
    summarise(
      presence_rate = mean(presence),
      n_days = n(),
      mean_detections_total = mean(n_detections_total, na.rm = TRUE),
      mean_species_detected = mean(n_species_detected, na.rm = TRUE),
      across(all_of(predictor_vars), first),
      .groups = "drop"
    ) %>%
    left_join(site_data %>% select(site_id, area_id, x, y), by = c("site_id", "area_id"))
}

apply_site_exclusion <- function(detection_data, site_data, site_exclusion) {
  if (!all(c("site_id", "remove") %in% names(site_exclusion))) {
    warning("site_exclusion.csv must contain site_id and remove columns. No sites were removed.")
    return(list(detection_data = detection_data, site_data = site_data))
  }
  
  exclusion_table <- site_exclusion %>%
    mutate(remove = as.logical(remove))
  
  sites_to_remove <- exclusion_table %>%
    filter(remove) %>%
    pull(site_id)
  
  list(
    detection_data = detection_data %>% filter(!site_id %in% sites_to_remove),
    site_data = site_data %>% filter(!site_id %in% sites_to_remove)
  )
}

assess_site_redundancy <- function(site_level_data, predictor_vars, site_exclusion) {
  if (!all(c("site_id", "remove") %in% names(site_exclusion))) {
    return(NULL)
  }
  
  dat <- site_level_data %>%
    drop_na(x, y)
  
  if (nrow(dat) < 2) return(NULL)
  
  coords <- dat %>%
    select(x, y) %>%
    as.matrix()
  
  spatial_dist <- as.matrix(dist(coords))
  rownames(spatial_dist) <- dat$site_id
  colnames(spatial_dist) <- dat$site_id
  diag(spatial_dist) <- NA
  
  env_scaled <- dat %>%
    select(all_of(predictor_vars)) %>%
    scale()
  
  environmental_dist <- as.matrix(dist(env_scaled))
  rownames(environmental_dist) <- dat$site_id
  colnames(environmental_dist) <- dat$site_id
  diag(environmental_dist) <- NA
  
  redundancy_table <- lapply(dat$site_id, function(id) {
    nearest_spatial <- names(which.min(spatial_dist[id, ]))
    nearest_environmental <- names(which.min(environmental_dist[id, ]))
    
    tibble(
      site_id = id,
      area_id = dat$area_id[dat$site_id == id],
      prevalence = dat$presence_rate[dat$site_id == id],
      mean_detections_total = dat$mean_detections_total[dat$site_id == id],
      mean_species_detected = dat$mean_species_detected[dat$site_id == id],
      nearest_spatial_site = nearest_spatial,
      spatial_distance = round(spatial_dist[id, nearest_spatial], 3),
      nearest_environmental_site = nearest_environmental,
      environmental_distance = round(environmental_dist[id, nearest_environmental], 3)
    )
  }) %>%
    bind_rows()
  
  redundancy_table %>%
    left_join(site_exclusion, by = "site_id") %>%
    mutate(
      remove = if_else(is.na(remove), FALSE, as.logical(remove)),
      decision = if_else(remove, "removed", "retained")
    ) %>%
    arrange(area_id, decision, site_id)
}

filter_predictors <- function(data, predictors, cutoff = 0.75) {
  x <- data %>%
    select(all_of(predictors)) %>%
    select(where(is.numeric))
  
  x <- x[, sapply(x, function(v) sd(v, na.rm = TRUE) > 0), drop = FALSE]
  
  if (ncol(x) < 2) {
    return(list(kept = names(x), removed = character(0), cor_mat = NULL))
  }
  
  cor_mat <- cor(x, use = "pairwise.complete.obs")
  remove_idx <- caret::findCorrelation(cor_mat, cutoff = cutoff)
  
  list(
    kept = names(x)[-remove_idx],
    removed = names(x)[remove_idx],
    cor_mat = cor_mat
  )
}

run_site_loocv_pls <- function(site_level_data, predictors) {
  dat <- site_level_data %>%
    select(site_id, area_id, presence_rate, x, y, all_of(predictors)) %>%
    drop_na()
  
  preds <- rep(NA_real_, nrow(dat))
  
  for (i in seq_len(nrow(dat))) {
    train_data <- dat[-i, ]
    test_data <- dat[i, ]
    ncomp_i <- min(3, length(predictors), nrow(train_data) - 1)
    
    if (ncomp_i < 1) next
    
    model <- train(
      presence_rate ~ .,
      data = train_data %>% select(-site_id, -area_id, -x, -y),
      method = "pls",
      trControl = trainControl(method = "none"),
      tuneGrid = data.frame(ncomp = ncomp_i)
    )
    
    preds[i] <- predict(model, newdata = test_data)
  }
  
  predictions <- dat %>%
    mutate(
      predicted = preds,
      residual = presence_rate - predicted
    )
  
  list(
    predictions = predictions,
    cor = cor(predictions$presence_rate, predictions$predicted, use = "complete.obs"),
    rmse = sqrt(mean((predictions$presence_rate - predictions$predicted)^2, na.rm = TRUE))
  )
}

run_leave_one_area_pls <- function(site_level_data, predictors) {
  dat <- site_level_data %>%
    select(site_id, area_id, presence_rate, all_of(predictors)) %>%
    drop_na()
  
  predictions <- map_dfr(unique(dat$area_id), function(area_name) {
    train_data <- dat %>% filter(area_id != area_name)
    test_data <- dat %>% filter(area_id == area_name)
    ncomp_i <- min(3, length(predictors), nrow(train_data) - 1)
    
    if (ncomp_i < 1) return(NULL)
    
    model <- train(
      presence_rate ~ .,
      data = train_data %>% select(-site_id, -area_id),
      method = "pls",
      trControl = trainControl(method = "none"),
      tuneGrid = data.frame(ncomp = ncomp_i)
    )
    
    test_data %>%
      mutate(
        predicted = predict(model, newdata = test_data),
        residual = presence_rate - predicted,
        left_out_area = area_name
      )
  })
  
  list(
    predictions = predictions,
    cor = cor(predictions$presence_rate, predictions$predicted, use = "complete.obs"),
    rmse = sqrt(mean((predictions$presence_rate - predictions$predicted)^2, na.rm = TRUE))
  )
}

run_moran_residuals <- function(prediction_data) {
  dat <- prediction_data %>%
    filter(!is.na(x), !is.na(y), !is.na(residual))
  
  coords <- as.matrix(dat[, c("x", "y")])
  distances <- as.matrix(dist(coords))
  diag(distances) <- NA
  
  weights <- 1 / distances
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

run_daily_pls_sdm <- function(data, predictors) {
  dat <- data %>%
    select(site_id, area_id, presence, all_of(predictors)) %>%
    drop_na()
  
  dat$presence <- factor(dat$presence, levels = c(0, 1), labels = c("absence", "presence"))
  
  if (length(unique(dat$presence)) < 2) {
    message("Skipping PLS: only one response class.")
    return(NULL)
  }
  
  folds <- groupKFold(group = dat$site_id, k = length(unique(dat$site_id)))
  
  ctrl <- trainControl(
    method = "cv",
    index = folds,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  
  tune_grid <- expand.grid(ncomp = 1:min(5, length(predictors)))
  
  tryCatch({
    model <- train(
      presence ~ .,
      data = dat %>% select(-site_id, -area_id),
      method = "pls",
      metric = "ROC",
      trControl = ctrl,
      tuneGrid = tune_grid
    )
    
    pred <- model$pred %>% filter(!is.na(presence))
    
    if (!"presence" %in% names(pred)) return(NULL)
    if (length(unique(pred$obs)) < 2) return(NULL)
    
    roc_obj <- pROC::roc(
      response = pred$obs,
      predictor = pred$presence,
      levels = c("absence", "presence"),
      quiet = TRUE
    )
    
    best_threshold <- pROC::coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = "threshold"
    )
    
    list(
      model = model,
      predictions = pred,
      auc = as.numeric(pROC::auc(roc_obj)),
      youden_threshold = as.numeric(best_threshold),
      predictors = predictors,
      data = dat
    )
  }, error = function(e) {
    message("PLS failed: ", e$message)
    return(NULL)
  })
}

run_area_specific_pls <- function(data, area_name, predictor_set, minimum_loggers = 3) {
  area_data <- data %>% filter(area_id == area_name)
  
  if (length(unique(area_data$site_id)) < minimum_loggers) {
    message("Skipping ", area_name, ": fewer than ", minimum_loggers, " sites.")
    return(NULL)
  }
  
  filtered <- filter_predictors(area_data, predictor_set, cutoff = correlation_cutoff)
  result <- run_daily_pls_sdm(area_data, filtered$kept)
  
  list(
    area_id = area_name,
    predictors = filtered$kept,
    removed_predictors = filtered$removed,
    model = result
  )
}

predict_sdm_raster_pls <- function(raster_stack, model_object, predictors) {
  missing_layers <- setdiff(predictors, names(raster_stack))
  
  if (length(missing_layers) > 0) {
    stop("Missing raster layers: ", paste(missing_layers, collapse = ", "))
  }
  
  predictor_raster <- raster_stack[[predictors]]
  names(predictor_raster) <- predictors
  
  prediction <- terra::predict(
    predictor_raster,
    model_object,
    fun = function(model, data) {
      data <- as.data.frame(data)
      probabilities <- predict(model, newdata = data, type = "prob")
      probabilities[, "presence"]
    },
    na.rm = TRUE
  )
  
  names(prediction) <- "suitability"
  prediction
}

make_extrapolation_map <- function(raster_stack, training_data, predictors) {
  predictor_raster <- raster_stack[[predictors]]
  flags <- vector("list", length(predictors))
  names(flags) <- predictors
  
  for (var in predictors) {
    train_min <- min(training_data[[var]], na.rm = TRUE)
    train_max <- max(training_data[[var]], na.rm = TRUE)
    flags[[var]] <- (predictor_raster[[var]] < train_min) | (predictor_raster[[var]] > train_max)
  }
  
  novelty <- Reduce(`+`, flags)
  names(novelty) <- "n_predictors_outside_training_range"
  novelty
}

prepared_detection_data <- prepare_detection_data(detection_data, breeding_start, breeding_end)
daily_data_all_sites <- build_daily_dataset(prepared_detection_data, target_species, predictor_vars)
site_level_all_sites <- build_site_dataset(daily_data_all_sites, site_data_raw, predictor_vars)

redundancy_table <- assess_site_redundancy(site_level_all_sites, predictor_vars, site_exclusion)

if (!is.null(redundancy_table)) {
  write.csv2(redundancy_table, file.path(table_dir, "Table_S1_site_redundancy_assessment.csv"), row.names = FALSE)
}

filtered_inputs <- apply_site_exclusion(prepared_detection_data, site_data_raw, site_exclusion)
prepared_detection_data <- filtered_inputs$detection_data
site_data <- filtered_inputs$site_data

daily_data <- build_daily_dataset(prepared_detection_data, target_species, predictor_vars)
site_level_data <- build_site_dataset(daily_data, site_data, predictor_vars)

period_datasets <- list(
  TOT = daily_data,
  PR = daily_data %>% filter(period == "PR"),
  NR = daily_data %>% filter(period == "NR")
)

main_data <- if (main_period == "TOT") {
  period_datasets$TOT
} else if (main_period == "PR") {
  period_datasets$PR
} else if (main_period == "NR") {
  period_datasets$NR
} else {
  stop("main_period must be one of: TOT, PR, NR")
}

main_site_level_data <- build_site_dataset(main_data, site_data, predictor_vars)

acoustic_summary <- daily_data %>%
  group_by(area_id, site_id) %>%
  summarise(
    n_days = n(),
    presences = sum(presence),
    absences = sum(presence == 0),
    prevalence = mean(presence),
    mean_detections_total = mean(n_detections_total, na.rm = TRUE),
    mean_species_detected = mean(n_species_detected, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(acoustic_summary, file.path(table_dir, "site_acoustic_sampling_summary.csv"), row.names = FALSE)

site_filter <- filter_predictors(main_site_level_data, predictor_vars, cutoff = correlation_cutoff)
site_loocv <- run_site_loocv_pls(main_site_level_data, site_filter$kept)
area_transfer <- run_leave_one_area_pls(main_site_level_data, site_filter$kept)
moran_residuals <- run_moran_residuals(site_loocv$predictions)

write.csv(site_loocv$predictions, file.path(table_dir, "site_level_LOOCV_predictions.csv"), row.names = FALSE)
write.csv(area_transfer$predictions, file.path(table_dir, "leave_one_area_out_predictions.csv"), row.names = FALSE)
write.csv(moran_residuals, file.path(table_dir, "moran_residuals_site_level_LOOCV.csv"), row.names = FALSE)

period_results <- map(period_datasets, function(period_data) {
  filtered <- filter_predictors(period_data, predictor_vars, cutoff = correlation_cutoff)
  run_daily_pls_sdm(period_data, filtered$kept)
})

period_auc_summary <- map_dfr(names(period_results), function(period_name) {
  result <- period_results[[period_name]]
  if (is.null(result)) return(NULL)
  
  tibble(
    period = period_name,
    model = "pls",
    auc = result$auc,
    youden_threshold = result$youden_threshold,
    n_predictors = length(result$predictors),
    predictors = paste(result$predictors, collapse = "; ")
  )
})

write.csv(period_auc_summary, file.path(table_dir, "period_level_AUC_summary.csv"), row.names = FALSE)

area_models <- map(
  names(area_rasters),
  ~ run_area_specific_pls(
    data = main_data,
    area_name = .x,
    predictor_set = predictor_vars,
    minimum_loggers = minimum_loggers_per_area
  )
)

names(area_models) <- names(area_rasters)

map_outputs <- list()
novelty_outputs <- list()

for (area_name in names(area_models)) {
  area_model <- area_models[[area_name]]
  
  if (is.null(area_model)) next
  if (is.null(area_model$model)) next
  
  suitability <- predict_sdm_raster_pls(
    raster_stack = area_rasters[[area_name]],
    model_object = area_model$model$model,
    predictors = area_model$model$predictors
  )
  
  suitability_file <- file.path(raster_dir, paste0(area_name, "_pls_suitability.tif"))
  writeRaster(suitability, suitability_file, overwrite = TRUE)
  map_outputs[[area_name]] <- suitability
  
  training_data <- main_data %>% filter(area_id == area_name)
  
  novelty <- make_extrapolation_map(
    raster_stack = area_rasters[[area_name]],
    training_data = training_data,
    predictors = area_model$model$predictors
  )
  
  novelty_file <- file.path(raster_dir, paste0(area_name, "_extrapolation_novelty.tif"))
  writeRaster(novelty, novelty_file, overwrite = TRUE)
  novelty_outputs[[area_name]] <- novelty
}

area_auc_summary <- map_dfr(names(area_models), function(area_name) {
  area_model <- area_models[[area_name]]
  
  if (is.null(area_model)) return(NULL)
  if (is.null(area_model$model)) return(NULL)
  
  tibble(
    area_id = area_name,
    period = main_period,
    model = "pls",
    auc = area_model$model$auc,
    youden_threshold = area_model$model$youden_threshold,
    n_sites = length(unique(area_model$model$data$site_id)),
    n_observations = nrow(area_model$model$data),
    predictors = paste(area_model$model$predictors, collapse = "; ")
  )
})

variable_importance <- map_dfr(names(area_models), function(area_name) {
  area_model <- area_models[[area_name]]
  
  if (is.null(area_model)) return(NULL)
  if (is.null(area_model$model)) return(NULL)
  
  importance <- varImp(area_model$model$model)$importance
  importance$predictor <- rownames(importance)
  rownames(importance) <- NULL
  
  importance %>%
    mutate(
      area_id = area_name,
      period = main_period,
      model = "pls"
    )
})

diagnostic_summary <- tibble(
  analysis = c(
    "site_level_leave_one_site_out_PLS",
    "site_level_leave_one_area_out_PLS",
    "moran_residuals_site_level_PLS"
  ),
  statistic = c("correlation", "correlation", "Moran_I"),
  value = c(site_loocv$cor, area_transfer$cor, moran_residuals$observed),
  rmse = c(site_loocv$rmse, area_transfer$rmse, NA_real_),
  p_value = c(NA_real_, NA_real_, moran_residuals$p_value)
)

paper_acoustic_summary <- acoustic_summary %>%
  mutate(
    prevalence = round(prevalence, 3),
    mean_detections_total = round(mean_detections_total, 1),
    mean_species_detected = round(mean_species_detected, 1)
  ) %>%
  arrange(area_id, site_id)

paper_period_auc <- period_auc_summary %>%
  mutate(
    auc = round(auc, 3),
    youden_threshold = round(youden_threshold, 3)
  )

paper_area_auc <- area_auc_summary %>%
  mutate(
    auc = round(auc, 3),
    youden_threshold = round(youden_threshold, 3)
  ) %>%
  arrange(area_id)

paper_diagnostics <- diagnostic_summary %>%
  mutate(
    value = round(value, 3),
    rmse = round(rmse, 3),
    p_value = signif(p_value, 3)
  )

write.csv(paper_acoustic_summary, file.path(table_dir, "Table_1_acoustic_sampling_summary.csv"), row.names = FALSE)
write.csv(paper_period_auc, file.path(table_dir, "Table_2_period_level_AUC_summary.csv"), row.names = FALSE)
write.csv(paper_area_auc, file.path(table_dir, "Table_3_area_specific_AUC_summary.csv"), row.names = FALSE)
write.csv(paper_diagnostics, file.path(table_dir, "Table_4_validation_diagnostics.csv"), row.names = FALSE)
write.csv(variable_importance, file.path(table_dir, "variable_importance_area_specific_PLS.csv"), row.names = FALSE)

write.csv2(paper_acoustic_summary, file.path(table_dir, "Table_1_acoustic_sampling_summary_excel.csv"), row.names = FALSE)
write.csv2(paper_period_auc, file.path(table_dir, "Table_2_period_level_AUC_summary_excel.csv"), row.names = FALSE)
write.csv2(paper_area_auc, file.path(table_dir, "Table_3_area_specific_AUC_summary_excel.csv"), row.names = FALSE)
write.csv2(paper_diagnostics, file.path(table_dir, "Table_4_validation_diagnostics_excel.csv"), row.names = FALSE)

p_site_loocv <- ggplot(site_loocv$predictions, aes(x = presence_rate, y = predicted, color = area_id)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    x = "Observed acoustic occurrence rate",
    y = "Predicted occurrence rate",
    color = "Area",
    title = "Site-level leave-one-site-out validation"
  ) +
  theme(
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

ggsave(file.path(figure_dir, "Figure_1_site_level_LOOCV.png"), p_site_loocv, width = 7, height = 5, dpi = 300)

p_area_transfer <- ggplot(area_transfer$predictions, aes(x = presence_rate, y = predicted, color = left_out_area)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    x = "Observed acoustic occurrence rate",
    y = "Predicted occurrence rate",
    color = "Left-out area",
    title = "Leave-one-area-out transferability test"
  ) +
  theme(
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

ggsave(file.path(figure_dir, "Figure_2_leave_one_area_out.png"), p_area_transfer, width = 7, height = 5, dpi = 300)

p_area_performance <- ggplot(paper_area_auc, aes(x = area_id, y = auc, fill = area_id)) +
  geom_col(width = 0.7) +
  theme_bw() +
  labs(
    x = "Area",
    y = "Leave-one-site-out AUC",
    fill = "Area",
    title = "Area-specific PLS model performance"
  ) +
  theme(
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

ggsave(file.path(figure_dir, "Figure_3_area_specific_model_performance.png"), p_area_performance, width = 8, height = 5, dpi = 300)

importance_plot_data <- variable_importance %>%
  rename(importance = Overall) %>%
  group_by(area_id) %>%
  slice_max(order_by = importance, n = 10, with_ties = FALSE) %>%
  ungroup()

p_importance <- ggplot(importance_plot_data, aes(x = reorder(predictor, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ area_id, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Predictor",
    y = "Variable importance",
    title = "Top predictors of acoustic suitability"
  ) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10)
  )

ggsave(file.path(figure_dir, "Figure_4_variable_importance_area_specific.png"), p_importance, width = 11, height = 6, dpi = 300)

if (length(map_outputs) > 0) {
  png(file.path(figure_dir, "Figure_5_suitability_maps.png"), width = 2200, height = 800, res = 200)
  par(mfrow = c(1, length(map_outputs)), mar = c(4, 4, 3, 5), oma = c(0, 0, 3, 0))
  walk(names(map_outputs), function(area_name) {
    plot(map_outputs[[area_name]], main = area_name, col = viridis(100), range = c(0, 1), cex.main = 1.2)
  })
  mtext("Predicted acoustic habitat suitability", outer = TRUE, cex = 1.2, font = 2, line = 1)
  dev.off()
}

if (length(novelty_outputs) > 0) {
  max_novelty <- max(map_dbl(novelty_outputs, ~ max(values(.x), na.rm = TRUE)), na.rm = TRUE)
  novelty_colours <- viridis(max_novelty + 1)
  
  png(file.path(figure_dir, "Figure_6_extrapolation_novelty_maps.png"), width = 2200, height = 800, res = 200)
  par(mfrow = c(1, length(novelty_outputs)), mar = c(4, 4, 3, 5), oma = c(0, 0, 3, 0))
  walk(names(novelty_outputs), function(area_name) {
    plot(
      novelty_outputs[[area_name]],
      main = area_name,
      col = novelty_colours,
      breaks = seq(-0.5, max_novelty + 0.5, by = 1),
      cex.main = 1.2
    )
  })
  mtext("Environmental novelty: predictors outside the training range", outer = TRUE, cex = 1.2, font = 2, line = 1)
  dev.off()
}

message("Pipeline completed.")
message("Outputs saved in: ", output_dir)
