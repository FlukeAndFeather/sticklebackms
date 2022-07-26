#!/usr/bin/env Rscript

library(furrr)
library(moments)
library(ranger)
library(RcppRoll)
library(rstickleback)
library(tidyverse)
library(zoo)

# Parameters --------------------------------------------------------------

test_local <- FALSE

# Rscript analysis/cluster/runstickleback.R 150 5 8 4 8 4 8 200 10 8

if (!test_local) {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args) != 10) {
    stop("Requires 10 arguments: win_size, tol, n_train, t_train, n_test, t_test, sb_trees, rf_trees, n_trials, n_cpu")
  }
  params <- list(
    win_size = as.integer(args[1]),
    tol = as.double(args[2]),
    n_train = as.integer(args[3]), # deployments
    t_train = as.double(args[4]), # hours
    n_test = as.integer(args[5]),
    t_test = as.double(args[6]),
    sb_trees = as.integer(args[7]),
    rf_trees = as.integer(args[8]),
    n_trials = as.integer(args[9]),
    n_cpu = as.integer(args[10])
  )
} else {
  params <- list(
    win_size = 150L,
    tol = 5,
    n_train = 8L, # deployments
    t_train = 4, # hours
    n_test = 8L,
    t_test = 4,
    sb_trees = 8L,
    rf_trees = 200L,
    n_trials = 10L,
    n_cpu = 2L
  )
}

# Read data ---------------------------------------------------------------

events <- arrow::read_parquet("analysis/data/raw_data/events.parquet") %>%
  filter(event == "lunge")
sensors <- arrow::read_parquet("analysis/data/raw_data/sensors.parquet")

# Fix an error where 5 deployments have missing speed values
missing_speed <- c("bw180828-49", "bw180830-42", "bw180830-48", "bw180905-42")
events <- filter(events, !deployid %in% missing_speed)
sensors <- filter(sensors, !deployid %in% missing_speed)

# Utility functions -------------------------------------------------------

create_features <- function(sensors, params) {
  win_size <- params$win_size
  roll_fn <- function(fn) {
    function(x) fn(x, n = win_size, fill = NA, align = "center")
  }
  roll_fn2 <- function(fn) {
    function(x) rollapply(x, width = win_size, FUN = fn, fill = NA, align = "center")
  }
  stat_fns <- list(
    min = roll_fn(roll_min),
    max = roll_fn(roll_max),
    mean = roll_fn(roll_mean),
    sd = roll_fn(roll_sd),
    skew = roll_fn2(skewness),
    kurt = roll_fn2(kurtosis)
  )
  features <- sensors %>%
    group_by(deployid) %>%
    mutate(across(c(depth, pitch, roll, speed),
                  stat_fns)) %>%
    ungroup() %>%
    left_join(events, by = c("deployid", "datetime")) %>%
    mutate(event = factor(ifelse(is.na(event), "non-event", "event"))) %>%
    drop_na()
}

durations <- function(sensors) {
  sensors %>%
    group_by(deployid) %>%
    summarize(d = max(datetime) - min(datetime)) %>%
    pull(d) %>%
    as.numeric(unit = "secs")
}

f1 <- function(tp, fp, fn) {
  tp / (tp + 1/2 * (fp + fn))
}

delta_r <- function(tp, fp, fn, d) {
  r = (tp + fn) / sum(d)
  r_hat = (tp + fp) / sum(d)
  abs((r - r_hat) / r)
}

assess_rf <- function(p, features, events, params) {
  tol <- params$tol
  features <- features %>%
    mutate(class = p) %>%
    filter(class == "event")
  assess_deployment <- function(f, e) {
    if (nrow(f) == 0) {
      e2 <- e %>%
        mutate(nearest = NA,
               error = Inf,
               outcome = "FN")
    } else {
      tryCatch({
        fdt <- c(min(f$datetime, e$datetime) - 1,
                 f$datetime,
                 max(f$datetime, e$datetime + 1))
        e2 <- e %>%
          mutate(nearest = fdt[findInterval(datetime, fdt)],
                 error = abs(as.numeric(datetime - nearest, unit = "secs")),
                 outcome = ifelse(error < tol, "TP", "FN"))
      }, error = function(err) browser())
    }
    f2 <- anti_join(f, e2, by = c("deployid", "datetime"))

    bind_rows(
      select(e2, deployid, datetime, outcome),
      transmute(f2, deployid, datetime, outcome = "FP")
    )
  }
  map_dfr(unique(events$deployid),
          ~ assess_deployment(filter(features, deployid == .x),
                              filter(events, deployid == .x)))
}

predict_rf <- function(rf_thr, newdat) {
  p <- predict(rf_thr$randomforest, newdat)
  ifelse(p$predictions[,1] >= rf_thr$thr, "event", "non-event")
}

## Split data into train/test ---------------------------------------------

split_data <- function(sensors, events, data_dir, i, params) {
  n_train <- params$n_train
  t_train <- params$t_train
  n_test <- params$n_test
  win_size <- params$win_size
  buffer <- win_size / 2 / 10 + 1

  # Randomly sample up to [t_train] hours of data from [n_train] deployments
  deployids <- sample(unique(events$deployid), size = n_train, replace = FALSE)
  get_start <- function(deployid) {
    e_t <- events$datetime[events$deployid == deployid]
    result <- sample(e_t, 1) - runif(1, buffer, t_train * 3600 - buffer)
  }
  start_times <- map(deployids, get_start)
  names(start_times) <- deployids
  train_sensors <- map_dfr(deployids, function(id) {
    start_time <- start_times[[id]]
    sensors %>%
      filter(deployid == id,
             datetime >= start_time,
             datetime <= start_time + 3600 * t_train)
  })
  set_nearest <- function(x, y) {
    # Assumes both are sorted
    stopifnot(!is.unsorted(x),
              !is.unsorted(y))
    # Assumes x within y
    stopifnot(x[1] > y[1],
              x[length(x)] < y[length(y)])

    y[findInterval(x, y)]
  }
  train_events <- map_dfr(deployids, function(id) {
    start_time <- start_times[[id]]
    valid_times <- train_sensors$datetime[train_sensors$deployid == id]
    max_time <- max(train_sensors$datetime[train_sensors$deployid == id])
    events %>%
        filter(deployid == id,
               # Filter out events near boundaries
               datetime > start_time + buffer,
               datetime < max_time - buffer) %>%
        mutate(datetime = set_nearest(datetime, valid_times))
  })

  # Needs refactoring SO BAD
  test_deployids <- sample(
    setdiff(unique(events$deployid), deployids),
    size = n_test,
    replace = FALSE
  )
  start_times <- map(test_deployids, get_start)
  names(start_times) <- test_deployids
  test_sensors <- map_dfr(test_deployids, function(id) {
    start_time <- start_times[[id]]
    sensors %>%
      filter(deployid == id,
             datetime >= start_time,
             datetime <= start_time + 3600 * t_train)
  })
  test_events <- map_dfr(test_deployids, function(id) {
    start_time <- start_times[[id]]
    valid_times <- test_sensors$datetime[test_sensors$deployid == id]
    max_time <- max(test_sensors$datetime[test_sensors$deployid == id])
    events %>%
      filter(deployid == id,
             # Filter out events near boundaries
             datetime > start_time + buffer,
             datetime < max_time - buffer) %>%
      mutate(datetime = set_nearest(datetime, valid_times))
  })

  # Check results
  tryCatch({
    stopifnot(
      setequal(unique(test_sensors$deployid), unique(test_events$deployid)),
      setequal(unique(train_sensors$deployid), unique(train_events$deployid))
    )
    for (id in deployids) {
      train_sensors_range <- range(filter(train_sensors, deployid == id)$datetime)
      train_events_range <- range(filter(train_events, deployid == id)$datetime)
      stopifnot(train_events_range[1] > train_sensors_range[1] + buffer,
                train_events_range[2] < train_sensors_range[2] - buffer)
    }
    for (id in test_deployids) {
      test_sensors_range <- range(filter(test_sensors, deployid == id)$datetime)
      test_events_range <- range(filter(test_events, deployid == id)$datetime)
      stopifnot(test_events_range[1] > test_sensors_range[1] + buffer,
                test_events_range[2] < test_sensors_range[2] - buffer)
    }
  }, error = function(e) browser())

  trial_dir <- file.path(data_dir, sprintf("%02.f", i))
  dir.create(trial_dir)
  saveRDS(train_sensors, file.path(trial_dir, "train_sensors.rds"))
  saveRDS(train_events, file.path(trial_dir, "train_events.rds"))
  saveRDS(test_sensors, file.path(trial_dir, "test_sensors.rds"))
  saveRDS(test_events, file.path(trial_dir, "test_events.rds"))
  trial_dir
}

## Train models -----------------------------------------------------------

train_stickleback <- function(trial_dir, params) {
  sb_trees <- params$sb_trees
  win_size <- params$win_size
  tol <- params$tol

  tryCatch({
    sensors <- readRDS(file.path(trial_dir, "train_sensors.rds")) %>%
      Sensors("deployid",
              "datetime",
              c("depth", "pitch", "roll", "speed"))
    events <- readRDS(file.path(trial_dir, "train_events.rds")) %>%
      Events("deployid",
             "datetime")
    tsc <- compose_tsc(module = "interval_based",
                       algorithm = "TimeSeriesForestClassifier",
                       params = list(n_estimators = sb_trees),
                       columns = columns(sensors))
    sb <- Stickleback(tsc,
                      win_size = win_size,
                      tol = tol,
                      nth = 5,
                      n_folds = 4)
    sb_fit(sb, sensors, events)
  }, error = function(e) browser())
  sb
}

train_randomforest <- function(trial_dir, params, strategy = c("proba", "sample")) {
  strategy = match.arg(strategy)
  rf_trees <- params$rf_trees

  sensors <- readRDS(file.path(trial_dir, "train_sensors.rds"))
  events <- readRDS(file.path(trial_dir, "train_events.rds"))
  features <- create_features(sensors, params)
  feat_cols <- colnames(features)[grepl(".*_.*", colnames(features))]

  if (strategy == "proba") {
    train_deployid <- sample(
      unique(events$deployid),
      size = max(round(length(unique(events$deployid)) * 0.8), 1)
    )
    feat_train <- filter(features, deployid %in% train_deployid)
    feat_valid <- filter(features, !deployid %in% train_deployid)

    rf_form <- as.formula(sprintf("event ~ %s",
                                  paste(feat_cols, collapse = "+")))
    m <- ranger(rf_form, feat_train, num.trees = rf_trees, probability = TRUE)

    pred <- predict(m, feat_valid)

    max_prob <- max(pred$predictions[,1])
    if (max_prob > 0) {
      thr <- seq(0, max_prob, by = 0.01)[-1]
      thr_f1 <- function(thr) {
        p <- ifelse(pred$predictions[, 1] >= thr, "event", "non-event")
        o <- assess_rf(p, feat_valid, events, params)
        f1(sum(o$outcome == "TP"),
           sum(o$outcome == "FP"),
           sum(o$outcome == "FN"))
      }
      valid_f1 <- map_dbl(thr, thr_f1)
      thr_hat <- thr[which.max(valid_f1)]
    } else {
      thr_hat <- 0.
    }
  } else if (strategy == "sample") {
    feat_events <- semi_join(features, events, by = c("deployid", "datetime"))
    feat_nonevents <- features %>%
      anti_join(feat_events, by = c("deployid", "datetime")) %>%
      sample_n(nrow(feat_events))
    feat_train <- rbind(feat_events, feat_nonevents)
    rf_form <- as.formula(sprintf("event ~ %s",
                                  paste(feat_cols, collapse = "+")))
    m <- ranger(rf_form, feat_train, num.trees = rf_trees, probability = TRUE)
    thr_hat <- 0.5
  }

  list(randomforest = m, thr = thr_hat)
}

test_stickleback <- function(m, trial_dir) {
  sensors <- readRDS(file.path(trial_dir, "test_sensors.rds"))
  events <- readRDS(file.path(trial_dir, "test_events.rds"))

  d <- durations(sensors)

  sensors <- Sensors(sensors,
                     "deployid",
                     "datetime",
                     c("depth", "pitch", "roll", "speed"))
  events <- Events(events,
                   "deployid",
                   "datetime")

  p <- sb_predict(m, sensors)
  o <- sb_assess(m, p, events) %>%
    as.data.frame()

  saveRDS(o, file.path(trial_dir, "_sbpredictions.rds"))

  list(f1 = f1(sum(o$outcome == "TP"),
               sum(o$outcome == "FP"),
               sum(o$outcome == "FN")),
       delta_r = delta_r(sum(o$outcome == "TP"),
                         sum(o$outcome == "FP"),
                         sum(o$outcome == "FN"),
                         d))
}

test_randomforest<- function(m, trial_dir, params) {
  sensors <- readRDS(file.path(trial_dir, "test_sensors.rds"))
  events <- readRDS(file.path(trial_dir, "test_events.rds"))

  d <- durations(sensors)
  features <- create_features(sensors, params)

  p <- predict_rf(m, features)
  o <- assess_rf(p, features, events, params)

  saveRDS(o, file.path(trial_dir, "_rfpredictions.rds"))

  list(f1 = f1(sum(o$outcome == "TP"),
               sum(o$outcome == "FP"),
               sum(o$outcome == "FN")),
       delta_r = delta_r(sum(o$outcome == "TP"),
                         sum(o$outcome == "FP"),
                         sum(o$outcome == "FN"),
                         d))
}

# Run cross validation ----------------------------------------------------

cv_trial <- function(i, sensors, events, data_dir, params) {
  # Split data
  trial_dir <- split_data(sensors, events, data_dir, i, params)

  # Train models
  sb <- train_stickleback(trial_dir, params)
  rf <- train_randomforest(trial_dir, params)

  # Test models
  sb_results <- test_stickleback(sb, trial_dir)
  rf_results <- test_randomforest(rf, trial_dir, params)

  # Results
  tibble(
    trial = i,
    sb_f1 = sb_results$f1,
    sb_delta_r = sb_results$delta_r,
    rf_f1 = rf_results$f1,
    rf_delta_r = rf_results$delta_r,
    rf_thr = rf$thr
  )
}

run_trials <- function(n_trials, parallel = FALSE) {
  data_dir <- sprintf("analysis/data/derived_data/cvtrials/%s",
                      format(Sys.time(), "%Y%m%d%H%M"))
  dir.create(data_dir)
  writeLines(sprintf("%s = %s", names(params), as.character(params)),
             file.path(data_dir, "params.txt"))
  map_fn <- if (parallel) future_map_dfr else (map_dfr)
  results <- map_fn(seq(params$n_trials),
                    cv_trial,
                    sensors = sensors,
                    events = events,
                    data_dir = data_dir,
                    params = params)
  write_csv(results, file.path(data_dir, "_results.csv"))
  results
}

print(paste("Start:", Sys.time()))
plan(multisession, workers = params$n_cpu)
run_trials(params$n_trials, parallel = TRUE)
print(paste("Finish:", Sys.time()))
