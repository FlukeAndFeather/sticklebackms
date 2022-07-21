#!/usr/bin/env Rscript

library(moments)
library(ranger)
library(RcppRoll)
library(rstickleback)
library(tidyverse)
library(zoo)

# Parameters --------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
  stop("Requires 8 arguments: win_size, tol, n_train, t_train, sb_trees, rf_trees, n_trials")
}
win_size <- as.integer(args[1])
tol <- as.double(args[2])
n_train <- as.integer(args[3]) # deployments
t_train <- as.double(args[4]) # hours
sb_trees <- as.integer(args[5])
rf_trees <- as.integer(args[6])
n_trials <- as.integer(args[7])

# Read data ---------------------------------------------------------------

keep <- c("bw170815-21", "bw180830-48", "bw180827-53", "bw180830-52b")
events <- arrow::read_parquet("analysis/data/raw_data/events.parquet") %>%
  filter(event == "lunge",
         deployid %in% keep)
sensors <- arrow::read_parquet("analysis/data/raw_data/sensors.parquet") %>%
  filter(deployid %in% keep)

# Utility functions -------------------------------------------------------

create_features <- function(sensors) {
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

assess_rf <- function(p, features, sensors, events) {
  features <- features %>%
    mutate(class = p$predictions) %>%
    filter(class == "event")
  assess_deployment <- function(f, e) {
    if (nrow(f) == 0) {
      e2 <- e %>%
        mutate(nearest = NA,
               error = Inf,
               outcome = "fn")
    } else {
      e2 <- e %>%
        mutate(nearest = f$datetime[findInterval(datetime, f$datetime)],
               error = abs(as.numeric(datetime - nearest, unit = "secs")),
               outcome = ifelse(error < tol, "tp", "fn"))
    }
    f2 <- anti_join(f, e2, by = c("deployid", "datetime"))

    tibble(tp = sum(e2$outcome == "tp"),
           fp = nrow(f2),
           fn = sum(e2$outcome == "fn"))
  }
  map_dfr(unique(events$deployid),
          ~ assess_deployment(filter(features, deployid == .x),
                              filter(events, deployid == .x))) %>%
    summarize(across(c(tp, fp, fn), sum))
}

## Split data into train/test ---------------------------------------------

split_data <- function(sensors, events) {
  # Randomly sample up to two hours of data from five deployments
  deployids <- sample(unique(events$deployid), size = n_train, replace = FALSE)
  get_start <- function(deployid) {
    e_t <- events$datetime[events$deployid == deployid]
    # Buffer according to window size such that at least one event falls within
    # the sample window
    buffer <- win_size / 2 / 10 + 1
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
    valid_times <- sensors$datetime[sensors$deployid == id]
    events %>%
        filter(deployid == id,
               # Filter out events near boundaries
               datetime > start_time + win_size / 2 / 10,
               datetime < start_time + 3600 * t_train - win_size / 2 / 10) %>%
        mutate(datetime = set_nearest(datetime, valid_times))
  })

  test_sensors <- anti_join(sensors, train_sensors, by = "deployid")
  test_events <- anti_join(events, train_events, by = "deployid") %>%
    group_by(deployid) %>%
    mutate(
      datetime = set_nearest(
        datetime,
        test_sensors$datetime[test_sensors$deployid == deployid[1]]
      )
    )
  list(train_sensors = train_sensors,
       train_events = train_events,
       test_sensors = test_sensors,
       test_events = test_events)
}

## Train models -----------------------------------------------------------

train_stickleback <- function(sensors, events) {
  sensors <- Sensors(sensors,
                     "deployid",
                     "datetime",
                     c("depth", "pitch", "roll", "speed"))
  events <- Events(events,
                   "deployid",
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

train_randomforest <- function(sensors, events) {
  features <- create_features(sensors)
  feat_cols <- colnames(features)[grepl(".*_.*", colnames(features))]
  rf_form <- as.formula(sprintf("event ~ %s",
                                paste(feat_cols, collapse = "+")))
  ranger(rf_form, features, num.trees = rf_trees)
}

test_stickleback <- function(m, sensors, events) {
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
    as.data.frame() %>%
    group_by(deployid) %>%
    summarize(tp = sum(outcome == "TP"),
              fp = sum(outcome == "FP"),
              fn = sum(outcome == "FN"))

  list(f1 = with(o, f1(tp, fp, fn)),
       delta_r = with(o, delta_r(tp, fp, fn, d)))
}

test_randomforest<- function(m, sensors, events) {
  d <- durations(sensors)
  features <- create_features(sensors)

  p <- predict(m, features)
  o <- assess_rf(p, features, sensors, events)

  list(f1 = with(o, f1(tp, fp, fn)),
       delta_r = with(o, delta_r(tp, fp, fn, d)))
}

# Run cross validation ----------------------------------------------------

run_cv <- function(i, sensors, events) {
  train_test <- split_data(sensors, events)

  # Train models
  sb <- train_stickleback(train_test$train_sensors, train_test$train_events)
  rf <- train_randomforest(train_test$train_sensors, train_test$train_events)

  # Make and assess predictions
  sb_results <- test_stickleback(sb,
                                 train_test$test_sensors,
                                 train_test$test_events)
  rf_results <- test_randomforest(rf,
                                  train_test$test_sensors,
                                  train_test$test_events)

  # Results
  tibble(
    trial = i,
    train_deployids = paste(unique(train_test$train_events$deployid),
                            collapse = ", "),
    sb_f1 = sb_results$f1,
    sb_delta_r = sb_results$delta_r,
    rf_f1 = rf_results$f1,
    rf_delta_r = rf_results$delta_r,
  )
}

print(paste("Start:", Sys.time()))
results <- map_dfr(seq(n_trials), run_cv, sensors = sensors, events = events)
print(paste("Finish:", Sys.time()))
write_csv(results, "analysis/data/derived_data/crossval.csv")