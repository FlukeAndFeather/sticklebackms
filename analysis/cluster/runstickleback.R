#!/usr/bin/env Rscript

# Slurm uses system Python, which causes an issue in `library(rstickleback)`
tryCatch(
  reticulate::conda_version(),
  error = function(e) reticulate::use_python(
    system2("which", "python3", stdout = TRUE)
  )
)

library(furrr)
library(moments)
library(ranger)
library(RcppRoll)
library(rstickleback)
library(tidyverse)
library(zoo)

# Parameters --------------------------------------------------------------

test_local <- TRUE

# Rscript analysis/cluster/runstickleback.R breath 150 5 8 4 8 4 8 200 10 8

if (!test_local) {
  args <- commandArgs(trailingOnly=TRUE)
} else {
  args <- c("lunge", "50", "5", "2", "1", "2", "1", "4", "100", "2", "1")
}
param_names <- c("behavior", "win_size", "tol", "n_train", "t_train",
                 "n_test", "t_test", "sb_trees", "rf_trees", "n_trials",
                 "n_cpu")

if (length(args) != length(param_names)) {
  stop(sprintf("Requires %i arguments: %s",
               length(param_names),
               paste(param_names, collapse = ", ")))
}

param_types <- str_split("cidididiiii", "")[[1]]
params <- map2(args, param_types, function(a, t) {
  switch (t,
          c = as.character(a),
          d = as.double(a),
          i = as.integer(a))
}) %>%
  set_names(param_names)

# Read data ---------------------------------------------------------------

events <- arrow::read_parquet("analysis/data/raw_data/events.parquet") %>%
  filter(event == params$behavior)
sensors <- arrow::read_parquet("analysis/data/raw_data/sensors.parquet")

# Fix an error where 4 deployments have missing speed values
missing_speed <- c("bw180828-49", "bw180830-42", "bw180830-48", "bw180905-42")
events <- filter(events, !deployid %in% missing_speed)
sensors <- filter(sensors, !deployid %in% missing_speed)

# Utility functions -------------------------------------------------------

create_features <- function(sensors, win_size) {
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
  feat_cols <- if (params$behavior == "lunge") {
    c("depth", "pitch", "roll", "speed")
  } else if (params$behavior == "breath") {
    c("depth", "pitch", "roll", "jerk")
  } else {
    stop(sprintf("Unknown behavior '%s'", params$behavior))
  }
  features <- sensors %>%
    group_by(deployid) %>%
    mutate(across(all_of(feat_cols), stat_fns)) %>%
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

fit_rf <- function(n_trees, win_size, sensors, events) {
  # If a window overlaps an event call it an event
  assign_class <- function(sensors_dt, events_dt, window) {
    overlaps <- logical(length(sensors_dt))
    sensors_secs <- as.numeric(sensors_dt - min(sensors_dt), unit = "secs")
    events_secs <- as.numeric(events_dt - min(sensors_dt), unit = "secs")
    w2 <- as.integer(window / 2)
    event_window <- function(e) {
      i <- as.integer(e * 10) + 1
      (i - w2):(i + w2)
    }
    event_idx <- map(events_secs, event_window) %>%
      unlist() %>%
      pmax(1) %>%
      pmin(length(overlaps))
    overlaps[event_idx] <- TRUE
    factor(ifelse(overlaps, "event", "non-event"),
           levels = c("event", "non-event"))
  }

  features <- create_features(sensors, win_size) %>%
    group_by(deployid) %>%
    mutate(class = assign_class(datetime,
                                events$datetime[events$deployid == deployid[1]],
                                win_size)) %>%
    ungroup()

  feat_train <- rbind(
    semi_join(features, events, by = c("deployid", "datetime")),
    sample_n(filter(features, class == "non-event"), size = nrow(events) * 10)
  )
  feat_cols <- colnames(features)[grepl(".*_.*", colnames(features))]

  rf_form <- as.formula(sprintf("event ~ %s",
                                paste(feat_cols, collapse = "+")))
  ranger(rf_form,
         feat_train,
         num.trees = n_trees)
}

optim_rf <- function(sensors, events, win_size, n_trees, tol) {
  deployids <- unique(events$deployid)
  num_test <- max(1, as.integer(length(deployids) * 0.3))
  test_deployid <- sample(deployids, num_test)
  train_deployid <- setdiff(deployids, test_deployid)

  train_sensors <- filter(sensors, !deployid %in% test_deployid)
  train_events <- filter(events, !deployid %in% test_deployid)
  test_sensors <- filter(sensors, deployid %in% test_deployid)
  test_events <- filter(events, deployid %in% test_deployid)

  rf_f1 <- function(win_size) {
    m <- fit_rf(n_trees, win_size, train_sensors, train_events)
    feat_test <- create_features(test_sensors, win_size)
    predictions <- predict(m, feat_test)
    o <- assess_rf(predictions, feat_test, test_events, tol)
    f1(sum(o$outcome == "TP"),
       sum(o$outcome == "FP"),
       sum(o$outcome == "FN"))
  }

  win_opt <- seq(log(win_size / 10),
                 log(win_size),
                 length.out = 10) %>%
    exp() %>%
    as.integer()

  par_f1 <- map_dbl(win_opt, rf_f1)

  win_opt[which.max(par_f1)]
}

assess_rf <- function(p, features, events, tol) {
  features <- features %>%
    mutate(class = p$predictions) %>%
    filter(class == "event")
  assess_deployment <- function(f, e) {
    if (nrow(f) == 0) {
      transmute(e, deployid, datetime, outcome = "FN")
    } else {
      fdt <- c(min(f$datetime, e$datetime) - 1,
               f$datetime,
               max(f$datetime, e$datetime + 1))
      e2 <- e %>%
        mutate(nearest = fdt[findInterval(datetime, fdt)],
               error = abs(as.numeric(datetime - nearest, unit = "secs")),
               outcome = ifelse(error < tol, "TP", "FN"))
      f2 <- anti_join(f, e2, by = c("deployid", "datetime"))
      bind_rows(
        select(e2, deployid, datetime, outcome),
        transmute(f2, deployid, datetime, outcome = "FP")
      )
    }
  }
  map_dfr(unique(events$deployid),
          ~ assess_deployment(filter(features, deployid == .x),
                              filter(events, deployid == .x)))
}

## Split data into train/test ---------------------------------------------

split_data <- function(sensors, events, data_dir, i, params) {
  n_train <- params$n_train
  t_train <- params$t_train
  n_test <- params$n_test
  win_size <- params$win_size
  buffer <- win_size / 2 / 10 + 1

  # Randomly sample up to [t_train] hours of data from [n_train] deployments
  train_ids <- sample(unique(events$deployid), size = n_train, replace = FALSE)
  get_start <- function(id) {
    e_dt <- sample(events$datetime[events$deployid == id], 1)
    offset <- runif(1, buffer, t_train * 3600 - buffer)
    e_dt - offset
  }
  train_start <- map(train_ids, get_start)
  names(train_start) <- train_ids
  train_sensors <- map_dfr(train_ids, function(id) {
    start_time <- train_start[[id]]
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
  train_events <- map_dfr(train_ids, function(id) {
    sensor_dt <- train_sensors$datetime[train_sensors$deployid == id]
    events %>%
        filter(deployid == id,
               # Filter out events near boundaries
               datetime > min(sensor_dt) + buffer,
               datetime < max(sensor_dt) - buffer) %>%
        mutate(datetime = set_nearest(datetime, sensor_dt))
  })

  # Needs refactoring SO BAD
  test_ids <- sample(
    setdiff(unique(events$deployid), train_ids),
    size = n_test,
    replace = FALSE
  )
  test_start <- map(test_ids, get_start)
  names(test_start) <- test_ids
  test_sensors <- map_dfr(test_ids, function(id) {
    start_time <- test_start[[id]]
    sensors %>%
      filter(deployid == id,
             datetime >= start_time,
             datetime <= start_time + 3600 * t_train)
  })
  test_events <- map_dfr(test_ids, function(id) {
    sensor_dt <- test_sensors$datetime[test_sensors$deployid == id]
    events %>%
      filter(deployid == id,
             # Filter out events near boundaries
             datetime > min(sensor_dt) + buffer,
             datetime < max(sensor_dt) - buffer) %>%
      mutate(datetime = set_nearest(datetime, sensor_dt))
  })

  # Check results
  tryCatch({
    stopifnot(
      setequal(unique(test_sensors$deployid), unique(test_events$deployid)),
      setequal(unique(train_sensors$deployid), unique(train_events$deployid))
    )
    for (id in train_ids) {
      train_sensors_range <- range(filter(train_sensors, deployid == id)$datetime)
      train_events_range <- range(filter(train_events, deployid == id)$datetime)
      stopifnot(train_events_range[1] > train_sensors_range[1] + buffer,
                train_events_range[2] < train_sensors_range[2] - buffer)
    }
    for (id in test_ids) {
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

  sb_vars <- if (params$behavior == "lunge") {
    c("depth", "pitch", "roll", "speed")
  } else if (params$behavior == "breath") {
    c("depth", "pitch", "roll", "jerk")
  } else {
    stop(sprintf("Unknown behavior '%s'", params$behavior))
  }

  sensors <- readRDS(file.path(trial_dir, "train_sensors.rds")) %>%
    Sensors("deployid",
            "datetime",
            sb_vars)
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
  sb
}

train_randomforest <- function(trial_dir, params) {
  rf_trees <- params$rf_trees
  win_size <- params$win_size
  tol <- params$tol

  sensors <- readRDS(file.path(trial_dir, "train_sensors.rds"))
  events <- readRDS(file.path(trial_dir, "train_events.rds"))

  trees_rng <- c(rf_trees / 2, rf_trees * 2)
  win_rng <- c(win_size / 10, win_size)
  rf_par <- optim_rf(sensors, events, win_size, rf_trees, tol)

  list(
    m = fit_rf(n_trees = rf_trees, win_size = rf_par, sensors, events),
    win_size = rf_par
  )
}

test_stickleback <- function(m, trial_dir, params) {
  sensors <- readRDS(file.path(trial_dir, "test_sensors.rds"))
  events <- readRDS(file.path(trial_dir, "test_events.rds"))

  d <- durations(sensors)

  sb_vars <- if (params$behavior == "lunge") {
    c("depth", "pitch", "roll", "speed")
  } else if (params$behavior == "breath") {
    c("depth", "pitch", "roll", "jerk")
  } else {
    stop(sprintf("Unknown behavior '%s'", params$behavior))
  }
  sensors <- Sensors(sensors,
                     "deployid",
                     "datetime",
                     sb_vars)
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
  features <- create_features(sensors, m$win_size)

  p <- predict(m$m, features)
  o <- assess_rf(p, features, events, params$tol)

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
  sb_results <- test_stickleback(sb, trial_dir, params)
  rf_results <- test_randomforest(rf, trial_dir, params)

  # Results
  tibble(
    trial = i,
    sb_f1 = sb_results$f1,
    sb_delta_r = sb_results$delta_r,
    rf_f1 = rf_results$f1,
    rf_delta_r = rf_results$delta_r,
    rf_winsz = rf$win_size
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
run_trials(params$n_trials, parallel = FALSE)
print(paste("Finish:", Sys.time()))
