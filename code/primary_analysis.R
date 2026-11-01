# Fit and evaluate CT prediction models and summarize missingness and domain shift.
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmnet)
  library(haven)
  library(lme4)
  library(mice)
  library(naniar)
  library(pROC)
  library(RANN)
  library(ranger)
  library(tidyr)
})

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1) args[[1]] else
  "data/raw"
output_dir <- if (length(args) >= 2) args[[2]] else "results_reproduced"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

seed_base <- 20260812L
pv_names <- paste0("PV", 1:5, "CT")
student_predictors <- c(
  "S_AGE", "S_SPECEFF", "S_GENEFF", "S_USECOM", "S_USEINF",
  "S_ACCONT", "S_SPECLASS", "S_GENCLASS", "S_CODLRN", "S_ICTFUT",
  "S_ICTLRN", "S_ICTNEG", "S_ICTPOS", "S_GENACT", "S_USESTD",
  "S_NISB", "PV1CIL"
)
school_predictors <- c(
  "P_RATTCH", "P_EXPLRN", "P_PRIORH", "P_PRIORS", "P_VWICT",
  "P_ICTUSE", "P_ICTCOM", "P_EXPTCH"
)
all_predictors <- c(student_predictors, school_predictors)
demographics <- c("S_SEX", "S_AGE", "S_HISEI", "S_HOMLIT", "S_IMMIG")

country_map <- c(
  DEU = "Germany", DNK = "Denmark", FIN = "Finland", FRA = "France",
  KOR = "Korea", LUX = "Luxembourg", PRT = "Portugal",
  USA = "United States"
)
europe_train_codes <- c("DEU", "FIN", "FRA", "LUX", "PRT")

safe_numeric <- function(x) as.numeric(zap_labels(x))

read_country <- function(code) {
  sf <- file.path(raw_dir, paste0("BSG", code, "I2.sav"))
  pf <- file.path(raw_dir, paste0("BCG", code, "I2.sav"))
  if (!file.exists(sf) || !file.exists(pf)) {
    stop("Missing official data file(s) for ", code, " in ", raw_dir)
  }
  s_keep <- c(
    "IDCNTRY", "IDSCHOOL", "IDSTUD", "TOTWGTS", demographics,
    student_predictors, pv_names, paste0("PV", 2:5, "CIL")
  )
  p_keep <- c("IDSCHOOL", school_predictors)
  s <- read_sav(sf, col_select = all_of(unique(s_keep)))
  p <- read_sav(pf, col_select = all_of(p_keep))
  stopifnot(!anyDuplicated(p$IDSCHOOL))
  d <- left_join(s, p, by = "IDSCHOOL")
  stopifnot(nrow(d) == nrow(s))
  for (v in setdiff(names(d), c("IDCNTRY", "IDSCHOOL", "IDSTUD"))) {
    d[[v]] <- safe_numeric(d[[v]])
  }
  d$COUNTRY <- unname(country_map[[code]])
  d$COUNTRY_CODE <- code
  d$SCHOOL_UID <- paste(code, d$IDSCHOOL, sep = "_")
  d$ROW_UID <- paste(code, d$IDSTUD, sep = "_")
  d
}

message("Reading and merging official ICILS 2018 student and school files...")
country_data <- lapply(names(country_map), read_country)
names(country_data) <- names(country_map)
all_data <- bind_rows(country_data)

cil_mode <- Sys.getenv("CIL_MODE", "matched")
stopifnot(cil_mode %in% c("matched", "fixed"))
cil_col_for <- function(pv) paste0("PV", sub("^PV([0-9])CT$", "\\1", pv), "CIL")
pred_matrix <- function(d, vars, pv) {
  m <- as.matrix(d[, vars, drop = FALSE])
  if (cil_mode == "matched" && "PV1CIL" %in% vars) {
    stopifnot(all(is.finite(d[[cil_col_for(pv)]])))
    m[, "PV1CIL"] <- d[[cil_col_for(pv)]]
  }
  m
}

write_csv_safe <- function(x, name) {
  write.csv(x, file.path(output_dir, name), row.names = FALSE, na = "")
}

metric_values <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  err <- predicted - observed
  mse <- mean(err^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(err))
  denom <- sum((observed - mean(observed))^2)
  r2 <- if (denom > 0) 1 - sum(err^2) / denom else NA_real_
  c(
    n = length(observed), MSE = mse, RMSE = rmse, MAE = mae, R2 = r2,
    mean_error = mean(err), target_mean = mean(observed),
    target_sd = sd(observed), RMSE_in_ICILS_SD = rmse / 100,
    RMSE_over_target_SD = rmse / sd(observed)
  )
}

make_group_foldid <- function(groups, k = 5L, seed = seed_base) {
  ug <- unique(groups)
  set.seed(seed)
  ug <- sample(ug, length(ug))
  group_fold <- rep(seq_len(min(k, length(ug))), length.out = length(ug))
  group_fold[match(groups, ug)]
}

make_group_holdout <- function(groups, prop = 0.20, seed = seed_base) {
  ug <- unique(groups)
  set.seed(seed)
  test_groups <- sample(ug, max(1L, round(length(ug) * prop)))
  groups %in% test_groups
}

scale_train_test <- function(x_train, x_test) {
  center <- colMeans(x_train)
  scale <- apply(x_train, 2, sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(
    train = sweep(sweep(x_train, 2, center, "-"), 2, scale, "/"),
    test = sweep(sweep(x_test, 2, center, "-"), 2, scale, "/"),
    center = center, scale = scale
  )
}

knn_predict <- function(x_train, y_train, x_new, k) {
  k <- min(k, nrow(x_train))
  idx <- RANN::nn2(data = x_train, query = x_new, k = k)$nn.idx
  if (is.null(dim(idx))) idx <- matrix(idx, ncol = 1)
  matrix_y <- matrix(y_train[idx], nrow = nrow(idx), ncol = ncol(idx))
  rowMeans(matrix_y)
}

fit_model <- function(model_name, x, y, groups, seed = seed_base) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  foldid <- make_group_foldid(groups, k = 5L, seed = seed)

  if (model_name %in% c("Lasso", "Ridge")) {
    alpha <- ifelse(model_name == "Lasso", 1, 0)
    set.seed(seed)
    cvfit <- cv.glmnet(
      x, y, alpha = alpha, foldid = foldid, standardize = TRUE,
      type.measure = "mse", nlambda = 60
    )
    return(list(
      model = model_name, object = cvfit$glmnet.fit,
      lambda = cvfit$lambda.min, hyperparameter = cvfit$lambda.min
    ))
  }

  if (model_name == "PCR") {
    max_comp <- min(ncol(x), nrow(x) - 2L)
    candidates <- seq_len(max_comp)
    fold_mse <- matrix(NA_real_, nrow = max(foldid), ncol = length(candidates))
    for (f in seq_len(max(foldid))) {
      tr <- foldid != f
      va <- !tr
      sc <- scale_train_test(x[tr, , drop = FALSE], x[va, , drop = FALSE])
      pc <- prcomp(sc$train, center = FALSE, scale. = FALSE, rank. = max_comp)
      tr_scores <- pc$x
      va_scores <- sc$test %*% pc$rotation
      for (j in seq_along(candidates)) {
        nc <- candidates[[j]]
        design <- cbind(1, tr_scores[, seq_len(nc), drop = FALSE])
        coef <- qr.solve(design, y[tr])
        pred <- cbind(1, va_scores[, seq_len(nc), drop = FALSE]) %*% coef
        fold_mse[f, j] <- mean((pred - y[va])^2)
      }
    }
    best_ncomp <- candidates[[which.min(colMeans(fold_mse))]]
    sc <- scale_train_test(x, x)
    pc <- prcomp(sc$train, center = FALSE, scale. = FALSE, rank. = best_ncomp)
    coef <- qr.solve(cbind(1, pc$x[, seq_len(best_ncomp), drop = FALSE]), y)
    return(list(
      model = model_name, object = pc, coef = coef,
      center = sc$center, scale = sc$scale, ncomp = best_ncomp,
      hyperparameter = best_ncomp
    ))
  }

  if (model_name == "Random forest") {
    candidates <- unique(pmax(1L, pmin(ncol(x), c(
      floor(sqrt(ncol(x))), floor(ncol(x) / 3), floor(ncol(x) / 2)
    ))))
    cv_mse <- rep(0, length(candidates))
    for (j in seq_along(candidates)) {
      fold_losses <- numeric(max(foldid))
      for (f in seq_len(max(foldid))) {
        tr <- foldid != f
        va <- !tr
        set.seed(seed + 100L * j + f)
        rf <- ranger(
          x = x[tr, , drop = FALSE], y = y[tr], mtry = candidates[[j]],
          num.trees = 200, min.node.size = 5, num.threads = 1,
          seed = seed + 100L * j + f
        )
        pred <- predict(rf, data = x[va, , drop = FALSE])$predictions
        fold_losses[[f]] <- mean((pred - y[va])^2)
      }
      cv_mse[[j]] <- mean(fold_losses)
    }
    best_mtry <- candidates[[which.min(cv_mse)]]
    set.seed(seed)
    rf <- ranger(
      x = x, y = y, mtry = best_mtry, num.trees = 500,
      min.node.size = 5, importance = "permutation", num.threads = 1,
      seed = seed
    )
    return(list(
      model = model_name, object = rf, hyperparameter = best_mtry
    ))
  }

  if (model_name == "KNN") {
    candidates <- c(5L, 11L, 21L, 31L)
    cv_mse <- rep(0, length(candidates))
    for (j in seq_along(candidates)) {
      fold_losses <- numeric(max(foldid))
      for (f in seq_len(max(foldid))) {
        tr <- foldid != f
        va <- !tr
        sc <- scale_train_test(x[tr, , drop = FALSE], x[va, , drop = FALSE])
        pred <- knn_predict(sc$train, y[tr], sc$test, candidates[[j]])
        fold_losses[[f]] <- mean((pred - y[va])^2)
      }
      cv_mse[[j]] <- mean(fold_losses)
    }
    best_k <- candidates[[which.min(cv_mse)]]
    sc <- scale_train_test(x, x)
    return(list(
      model = model_name, x = sc$train, y = y, center = sc$center,
      scale = sc$scale, k = best_k, hyperparameter = best_k
    ))
  }

  stop("Unsupported model: ", model_name)
}

predict_model <- function(fit, newx) {
  newx <- as.matrix(newx)
  if (fit$model %in% c("Lasso", "Ridge")) {
    return(as.numeric(predict(fit$object, newx = newx, s = fit$lambda)))
  }
  if (fit$model == "PCR") {
    sx <- sweep(sweep(newx, 2, fit$center, "-"), 2, fit$scale, "/")
    scores <- sx %*% fit$object$rotation[, seq_len(fit$ncomp), drop = FALSE]
    return(as.numeric(cbind(1, scores) %*% fit$coef))
  }
  if (fit$model == "Random forest") {
    return(as.numeric(predict(fit$object, data = newx)$predictions))
  }
  if (fit$model == "KNN") {
    sx <- sweep(sweep(newx, 2, fit$center, "-"), 2, fit$scale, "/")
    return(knn_predict(fit$x, fit$y, sx, fit$k))
  }
  stop("Unsupported fit")
}

permutation_importance <- function(fit, x, y, repeats = 3L, seed = seed_base, groups = NULL) {
  baseline <- mean((predict_model(fit, x) - y)^2)
  out <- vector("list", ncol(x))
  for (j in seq_len(ncol(x))) {
    losses <- numeric(repeats)
    for (r in seq_len(repeats)) {
      set.seed(seed + 1000L * j + r)
      xp <- x
      if (colnames(x)[j] %in% school_predictors && !is.null(groups)) {
        ug <- unique(groups)
        values <- x[match(ug, groups), j]
        xp[, j] <- sample(values)[match(groups, ug)]
      } else {
        xp[, j] <- sample(xp[, j])
      }
      losses[[r]] <- mean((predict_model(fit, xp) - y)^2)
    }
    out[[j]] <- data.frame(
      variable = colnames(x)[[j]], baseline_MSE = baseline,
      permuted_MSE = mean(losses), delta_MSE = mean(losses) - baseline,
      percent_MSE_increase = 100 * (mean(losses) - baseline) / baseline
    )
  }
  bind_rows(out)
}

models <- c("Lasso", "Ridge", "PCR", "Random forest", "KNN")

message("Computing missingness, sample-flow, demographic, MCAR, and ICC diagnostics...")
missing_rates <- bind_rows(lapply(names(country_data), function(code) {
  d <- country_data[[code]]
  bind_rows(lapply(c(all_predictors, pv_names), function(v) {
    data.frame(
      country = unname(country_map[[code]]), country_code = code,
      variable = v, n_total = nrow(d), n_missing = sum(is.na(d[[v]])),
      percent_missing = 100 * mean(is.na(d[[v]]))
    )
  }))
}))
write_csv_safe(missing_rates, "missing_rates_by_country_variable.csv")

sample_flow <- bind_rows(lapply(names(country_data), function(code) {
  d <- country_data[[code]]
  student_complete <- complete.cases(d[, c(student_predictors, pv_names)])
  school_complete <- complete.cases(d[, school_predictors])
  analytic <- complete.cases(d[, c(all_predictors, pv_names)])
  data.frame(
    country = unname(country_map[[code]]), country_code = code,
    raw_students = nrow(d), schools = n_distinct(d$SCHOOL_UID),
    all_CT_PVs_observed = sum(complete.cases(d[, pv_names])),
    student_side_complete = sum(student_complete),
    school_side_complete = sum(school_complete),
    complete_case_analytic_n = sum(analytic),
    retained_percent = 100 * mean(analytic)
  )
}))
write_csv_safe(sample_flow, "sample_flow_by_country.csv")

smd_continuous <- function(x, g) {
  x1 <- x[g & !is.na(x)]
  x0 <- x[!g & !is.na(x)]
  pooled <- sqrt((var(x1) + var(x0)) / 2)
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  (mean(x1) - mean(x0)) / pooled
}

retention_demographics <- bind_rows(lapply(names(country_data), function(code) {
  d <- country_data[[code]]
  retained <- complete.cases(d[, c(all_predictors, pv_names)])
  bind_rows(lapply(demographics, function(v) {
    data.frame(
      country = unname(country_map[[code]]), country_code = code,
      variable = v,
      retained_n_observed = sum(retained & !is.na(d[[v]])),
      excluded_n_observed = sum(!retained & !is.na(d[[v]])),
      retained_mean = mean(d[[v]][retained], na.rm = TRUE),
      excluded_mean = mean(d[[v]][!retained], na.rm = TRUE),
      standardized_mean_difference = smd_continuous(d[[v]], retained)
    )
  }))
}))
write_csv_safe(retention_demographics, "retained_vs_excluded_demographics.csv")

mcar_results <- bind_rows(lapply(names(country_data), function(code) {
  d <- country_data[[code]][, c(all_predictors, "PV1CT")]
  result <- tryCatch(naniar::mcar_test(d), error = function(e) e)
  if (inherits(result, "error")) {
    data.frame(
      country = unname(country_map[[code]]), country_code = code,
      statistic = NA_real_, df = NA_real_, p_value = NA_real_,
      status = paste("Not estimable:", conditionMessage(result))
    )
  } else {
    data.frame(
      country = unname(country_map[[code]]), country_code = code,
      statistic = result$statistic, df = result$df, p_value = result$p.value,
      status = "Estimated (Little's test of MCAR; does not test MAR)"
    )
  }
}))
write_csv_safe(mcar_results, "little_mcar_test_by_country.csv")

icc_results <- bind_rows(lapply(names(country_data), function(code) {
  d <- country_data[[code]]
  bind_rows(lapply(pv_names, function(pv) {
    dd <- d[complete.cases(d[, c("SCHOOL_UID", pv)]), ]
    fit <- lmer(as.formula(paste(pv, "~ 1 + (1 | SCHOOL_UID)")), data = dd, REML = TRUE)
    vc <- as.data.frame(VarCorr(fit))
    between <- vc$vcov[vc$grp == "SCHOOL_UID"][[1]]
    within <- vc$vcov[vc$grp == "Residual"][[1]]
    data.frame(
      country = unname(country_map[[code]]), country_code = code,
      plausible_value = pv, schools = n_distinct(dd$SCHOOL_UID), n = nrow(dd),
      between_school_variance = between, within_school_variance = within,
      ICC = between / (between + within)
    )
  }))
}))
write_csv_safe(icc_results, "school_icc_by_country_and_plausible_value.csv")

message("Running within-country school-blocked models across five plausible values...")
target_splits <- list()
within_metrics <- list()
within_predictions <- list()
importance_results <- list()
hyperparameters <- list()
counter <- 0L

for (code in c("KOR", "DNK")) {
  d <- country_data[[code]]
  keep <- complete.cases(d[, c(all_predictors, pv_names)])
  d <- d[keep, ]
  test_flag <- make_group_holdout(d$SCHOOL_UID, prop = 0.20, seed = seed_base + match(code, c("KOR", "DNK")))
  target_splits[[code]] <- list(data = d, test_flag = test_flag)

  for (feature_set in c("All predictors", "Excluding PV1CIL")) {
    predictors <- if (feature_set == "All predictors") all_predictors else setdiff(all_predictors, "PV1CIL")
    groups_train <- d$SCHOOL_UID[!test_flag]

    for (pv in pv_names) {
      x_train <- pred_matrix(d[!test_flag, ], predictors, pv)
      x_test <- pred_matrix(d[test_flag, ], predictors, pv)
      y_train <- d[[pv]][!test_flag]
      y_test <- d[[pv]][test_flag]
      for (model_name in models) {
        counter <- counter + 1L
        message("  within ", code, " | ", feature_set, " | ", pv, " | ", model_name)
        fit <- fit_model(model_name, x_train, y_train, groups_train, seed_base + counter)
        pred <- predict_model(fit, x_test)
        met <- as.data.frame(as.list(metric_values(y_test, pred)))
        within_metrics[[counter]] <- cbind(
          data.frame(
            context = "Within-country school-blocked holdout",
            target_country = unname(country_map[[code]]), target_code = code,
            feature_set = feature_set, plausible_value = pv,
            model = model_name, train_students = nrow(x_train),
            test_students = nrow(x_test), train_schools = n_distinct(groups_train),
            test_schools = n_distinct(d$SCHOOL_UID[test_flag]),
            hyperparameter = fit$hyperparameter
          ), met
        )
        hyperparameters[[length(hyperparameters) + 1L]] <- data.frame(
          context = "Within-country", target_code = code,
          feature_set = feature_set, plausible_value = pv, model = model_name,
          hyperparameter = fit$hyperparameter
        )
        within_predictions[[counter]] <- data.frame(
          context = "Within-country", target_country = unname(country_map[[code]]),
          target_code = code, feature_set = feature_set, plausible_value = pv,
          model = model_name, row_uid = d$ROW_UID[test_flag],
          school_uid = d$SCHOOL_UID[test_flag], observed = y_test,
          predicted = pred, residual = y_test - pred
        )
        if (feature_set == "All predictors" && pv == "PV1CT") {
          imp <- permutation_importance(
            fit, x_test, y_test, repeats = 3L,
            seed = seed_base + 5000L + counter, groups = d$SCHOOL_UID[test_flag]
          )
          imp$context <- "Within-country"
          imp$target_country <- unname(country_map[[code]])
          imp$target_code <- code
          imp$model <- model_name
          importance_results[[length(importance_results) + 1L]] <- imp
        }
      }
    }
  }
}

within_metrics_df <- bind_rows(within_metrics)
within_predictions_df <- bind_rows(within_predictions)
write_csv_safe(within_metrics_df, "within_country_metrics_all_plausible_values.csv")
write_csv_safe(within_predictions_df, "within_country_predictions_all_plausible_values.csv")

message("Running five-country European external-validation models...")
europe <- bind_rows(country_data[europe_train_codes])
europe <- europe[complete.cases(europe[, c(all_predictors, pv_names)]), ]
external_metrics <- list()
external_predictions <- list()
external_counter <- 0L

for (feature_set in c("All predictors", "Excluding PV1CIL")) {
  predictors <- if (feature_set == "All predictors") all_predictors else setdiff(all_predictors, "PV1CIL")
  groups_train <- europe$SCHOOL_UID
  for (pv in pv_names) {
    x_train <- pred_matrix(europe, predictors, pv)
    y_train <- europe[[pv]]
    for (model_name in models) {
      external_counter <- external_counter + 1L
      message("  external Europe | ", feature_set, " | ", pv, " | ", model_name)
      fit <- fit_model(
        model_name, x_train, y_train, groups_train,
        seed_base + 10000L + external_counter
      )
      hyperparameters[[length(hyperparameters) + 1L]] <- data.frame(
        context = "Europe-trained external", target_code = "KOR/DNK",
        feature_set = feature_set, plausible_value = pv, model = model_name,
        hyperparameter = fit$hyperparameter
      )
      for (code in c("KOR", "DNK")) {
        split <- target_splits[[code]]
        dtest <- split$data[split$test_flag, ]
        x_test <- pred_matrix(dtest, predictors, pv)
        y_test <- dtest[[pv]]
        pred <- predict_model(fit, x_test)
        met <- as.data.frame(as.list(metric_values(y_test, pred)))
        external_metrics[[length(external_metrics) + 1L]] <- cbind(
          data.frame(
            context = "Five-country Europe-trained external validation",
            target_country = unname(country_map[[code]]), target_code = code,
            feature_set = feature_set, plausible_value = pv,
            model = model_name, train_students = nrow(x_train),
            test_students = nrow(x_test), train_schools = n_distinct(groups_train),
            test_schools = n_distinct(dtest$SCHOOL_UID),
            hyperparameter = fit$hyperparameter
          ), met
        )
        external_predictions[[length(external_predictions) + 1L]] <- data.frame(
          context = "Europe-trained external", target_country = unname(country_map[[code]]),
          target_code = code, feature_set = feature_set, plausible_value = pv,
          model = model_name, row_uid = dtest$ROW_UID,
          school_uid = dtest$SCHOOL_UID, observed = y_test,
          predicted = pred, residual = y_test - pred
        )
        if (feature_set == "All predictors" && pv == "PV1CT") {
          imp <- permutation_importance(
            fit, x_test, y_test, repeats = 3L,
            seed = seed_base + 15000L + external_counter + match(code, c("KOR", "DNK")), groups = dtest$SCHOOL_UID
          )
          imp$context <- "Europe-trained external"
          imp$target_country <- unname(country_map[[code]])
          imp$target_code <- code
          imp$model <- model_name
          importance_results[[length(importance_results) + 1L]] <- imp
        }
      }
    }
  }
}

external_metrics_df <- bind_rows(external_metrics)
external_predictions_df <- bind_rows(external_predictions)
write_csv_safe(external_metrics_df, "external_validation_metrics_all_plausible_values.csv")
write_csv_safe(external_predictions_df, "external_validation_predictions_all_plausible_values.csv")

all_metrics <- bind_rows(within_metrics_df, external_metrics_df)
pooled_metrics <- all_metrics |>
  group_by(context, target_country, target_code, feature_set, model) |>
  summarise(
    plausible_values = n(), n = mean(n),
    MSE_mean = mean(MSE), MSE_SD_across_PVs = sd(MSE),
    RMSE_mean = mean(RMSE), RMSE_SD_across_PVs = sd(RMSE),
    MAE_mean = mean(MAE), MAE_SD_across_PVs = sd(MAE),
    R2_mean = mean(R2), R2_SD_across_PVs = sd(R2),
    mean_error_mean = mean(mean_error),
    RMSE_in_ICILS_SD_mean = mean(RMSE_in_ICILS_SD),
    RMSE_over_target_SD_mean = mean(RMSE_over_target_SD),
    .groups = "drop"
  )
write_csv_safe(pooled_metrics, "model_metrics_pooled_across_plausible_values.csv")

performance_delta <- pooled_metrics |>
  select(context, target_country, target_code, feature_set, model, RMSE_mean, R2_mean) |>
  pivot_wider(names_from = context, values_from = c(RMSE_mean, R2_mean)) |>
  mutate(
    external_minus_within_RMSE =
      `RMSE_mean_Five-country Europe-trained external validation` -
      `RMSE_mean_Within-country school-blocked holdout`,
    external_minus_within_R2 =
      `R2_mean_Five-country Europe-trained external validation` -
      `R2_mean_Within-country school-blocked holdout`
  )
write_csv_safe(performance_delta, "external_vs_within_performance_delta.csv")

importance_df <- bind_rows(importance_results) |>
  group_by(context, target_country, target_code, model) |>
  mutate(
    positive_delta = pmax(delta_MSE, 0),
    normalized_importance = ifelse(
      sum(positive_delta) > 0, positive_delta / sum(positive_delta), 0
    ),
    rank = rank(-delta_MSE, ties.method = "min")
  ) |>
  ungroup()
write_csv_safe(importance_df, "permutation_importance_all_models.csv")
write_csv_safe(bind_rows(hyperparameters), "selected_hyperparameters.csv")

domain_shift_pair <- function(train, target, predictors = all_predictors) {
  smd <- sapply(predictors, function(v) {
    pooled <- sqrt((var(train[[v]]) + var(target[[v]])) / 2)
    if (!is.finite(pooled) || pooled == 0) return(NA_real_)
    (mean(target[[v]]) - mean(train[[v]])) / pooled
  })
  cov_train <- cov(train[, predictors])
  inv_cov <- tryCatch(solve(cov_train), error = function(e) MASS::ginv(cov_train))
  mean_diff <- colMeans(target[, predictors]) - colMeans(train[, predictors])
  mahal <- as.numeric(sqrt(t(mean_diff) %*% inv_cov %*% mean_diff))

  n_bal <- min(nrow(train), nrow(target), 3000L)
  set.seed(seed_base)
  a <- train[sample(seq_len(nrow(train)), n_bal), predictors]
  b <- target[sample(seq_len(nrow(target)), n_bal), predictors]
  xx <- as.matrix(bind_rows(a, b))
  yy <- c(rep(0, n_bal), rep(1, n_bal))

  set.seed(seed_base)
  ia <- sample(seq_len(nrow(train)), n_bal)
  ib <- sample(seq_len(nrow(target)), n_bal)
  groups <- c(train$SCHOOL_UID[ia], target$SCHOOL_UID[ib])
  folds <- c(make_group_foldid(groups[yy == 0], seed=seed_base+1L),
             make_group_foldid(groups[yy == 1], seed=seed_base+2L))
  aucs <- numeric(5)
  for (f in 1:5) {
    fit <- cv.glmnet(
      xx[folds != f, ], yy[folds != f], family = "binomial", alpha = 0,
      foldid = make_group_foldid(groups[folds != f], seed=seed_base+f), standardize = TRUE
    )
    pr <- as.numeric(predict(fit, newx = xx[folds == f, ], s = "lambda.min", type = "response"))
    aucs[[f]] <- as.numeric(pROC::auc(yy[folds == f], pr, quiet = TRUE))
  }
  list(
    summary = data.frame(
      mean_absolute_SMD = mean(abs(smd), na.rm = TRUE),
      max_absolute_SMD = max(abs(smd), na.rm = TRUE),
      mahalanobis_mean_distance = mahal,
      domain_classifier_AUC = mean(aucs)
    ),
    variable = data.frame(variable = names(smd), standardized_mean_difference = smd)
  )
}

domain_summaries <- list()
domain_variables <- list()
for (code in c("KOR", "DNK")) {
  target <- target_splits[[code]]$data[target_splits[[code]]$test_flag, ]
  shift <- domain_shift_pair(europe, target)
  shift$summary$training_context <- "Five European countries"
  shift$summary$target_country <- unname(country_map[[code]])
  shift$summary$target_code <- code
  shift$summary$PV1CT_mean_difference <- mean(target$PV1CT) - mean(europe$PV1CT)
  shift$summary$PV1CT_SD_ratio <- sd(target$PV1CT) / sd(europe$PV1CT)
  shift$variable$training_context <- "Five European countries"
  shift$variable$target_country <- unname(country_map[[code]])
  shift$variable$target_code <- code
  domain_summaries[[code]] <- shift$summary
  domain_variables[[code]] <- shift$variable
}
write_csv_safe(bind_rows(domain_summaries), "domain_shift_summary_europe_to_targets.csv")
write_csv_safe(bind_rows(domain_variables), "domain_shift_by_variable_europe_to_targets.csv")

message("Running leave-one-country-out diagnostics for lasso and random forest...")
loco_metrics <- list()
loco_shift <- list()
loco_counter <- 0L
for (target_code in names(country_map)) {
  target <- country_data[[target_code]]
  target <- target[complete.cases(target[, c(all_predictors, "PV1CT")]), ]
  train_codes <- setdiff(names(country_map), target_code)
  train <- bind_rows(country_data[train_codes])
  train <- train[complete.cases(train[, c(all_predictors, "PV1CT")]), ]
  shift <- domain_shift_pair(train, target)
  shift$summary$target_country <- unname(country_map[[target_code]])
  shift$summary$target_code <- target_code
  loco_shift[[target_code]] <- shift$summary
  for (model_name in c("Lasso", "Random forest")) {
    loco_counter <- loco_counter + 1L
    message("  LOCO ", target_code, " | ", model_name)
    fit <- fit_model(
      model_name, as.matrix(train[, all_predictors]), train$PV1CT,
      train$SCHOOL_UID, seed_base + 20000L + loco_counter
    )
    pred <- predict_model(fit, as.matrix(target[, all_predictors]))
    met <- as.data.frame(as.list(metric_values(target$PV1CT, pred)))
    loco_metrics[[loco_counter]] <- cbind(
      data.frame(
        target_country = unname(country_map[[target_code]]),
        target_code = target_code, model = model_name,
        training_countries = paste(unname(country_map[train_codes]), collapse = "; "),
        train_n = nrow(train), target_n = nrow(target),
        hyperparameter = fit$hyperparameter
      ), met
    )
  }
}
loco_metrics_df <- bind_rows(loco_metrics)
loco_shift_df <- bind_rows(loco_shift)
write_csv_safe(loco_metrics_df, "leave_one_country_out_metrics_PV1CT.csv")
write_csv_safe(loco_shift_df, "leave_one_country_out_domain_shift.csv")

loco_association <- left_join(loco_metrics_df, loco_shift_df, by = c("target_country", "target_code")) |>
  group_by(model) |>
  summarise(
    countries = n(),
    spearman_RMSE_vs_mean_abs_SMD = cor(RMSE, mean_absolute_SMD, method = "spearman"),
    pearson_RMSE_vs_mean_abs_SMD = cor(RMSE, mean_absolute_SMD),
    spearman_RMSE_vs_domain_AUC = cor(RMSE, domain_classifier_AUC, method = "spearman"),
    .groups = "drop"
  )
write_csv_safe(loco_association, "domain_shift_performance_association.csv")

plot_perf <- pooled_metrics |>
  filter(feature_set == "All predictors") |>
  ggplot(aes(x = model, y = RMSE_mean, color = context, group = context)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.25)) +
  geom_line(position = position_dodge(width = 0.25)) +
  facet_wrap(~ target_country) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    x = NULL, y = "RMSE (ICILS CT scale points)", color = NULL,
    title = "Within-country and external-validation error",
    subtitle = "Means across five plausible values; school-blocked target holdouts"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(file.path(output_dir, "figures", "performance_comparison.png"), plot_perf, width = 8.5, height = 5.2, dpi = 300)

top_imp <- importance_df |>
  group_by(context, target_country, model) |>
  slice_max(delta_MSE, n = 8, with_ties = FALSE) |>
  ungroup() |>
  mutate(variable = reorder(variable, delta_MSE))
plot_imp <- ggplot(top_imp, aes(x = delta_MSE, y = variable, fill = model)) +
  geom_col(show.legend = FALSE) +
  facet_grid(target_country ~ context + model, scales = "free_y") +
  labs(
    x = "Increase in MSE after permutation", y = NULL,
    title = "Permutation importance by model and validation context",
    subtitle = "PV1CT; positive values indicate predictive contribution"
  ) +
  theme_minimal(base_size = 8)
ggsave(file.path(output_dir, "figures", "feature_importance.png"), plot_imp, width = 15, height = 8, dpi = 300)

plot_missing <- missing_rates |>
  mutate(variable = factor(variable, levels = rev(unique(variable)))) |>
  ggplot(aes(x = country, y = variable, fill = percent_missing)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  labs(
    x = NULL, y = NULL, fill = "% missing",
    title = "Missingness in candidate predictors and CT plausible values"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(output_dir, "figures", "missingness_heatmap.png"), plot_missing, width = 9, height = 8, dpi = 300)

variable_dictionary <- data.frame(
  variable = c(all_predictors, pv_names),
  level = c(rep("Student", length(student_predictors)), rep("School", length(school_predictors)), rep("Student outcome", length(pv_names))),
  role = c(
    rep("Predictor", length(student_predictors) - 1), "Proximal achievement covariate / sensitivity target",
    rep("Predictor", length(school_predictors)), rep("Computational thinking plausible value", length(pv_names))
  )
)

construct_map <- c(
  S_AGE = "Demographic background", S_NISB = "Socioeconomic background",
  S_SPECEFF = "Specialist ICT self-efficacy", S_GENEFF = "General ICT self-efficacy",
  S_USECOM = "ICT communication use", S_USEINF = "ICT information exchange",
  S_ACCONT = "Online content access", S_SPECLASS = "Specialist classroom ICT use",
  S_GENCLASS = "General classroom ICT use", S_CODLRN = "Coding learning opportunities",
  S_ICTFUT = "Future ICT expectations", S_ICTLRN = "ICT learning opportunities",
  S_ICTNEG = "Negative ICT attitudes", S_ICTPOS = "Positive ICT attitudes",
  S_GENACT = "General application use", S_USESTD = "ICT use for study",
  PV1CIL = "Computer and information literacy achievement",
  P_RATTCH = "School staffing/resources", P_EXPLRN = "Teacher ICT-use expectations",
  P_PRIORH = "Hardware/resource priorities", P_PRIORS = "Support/professional-learning priorities",
  P_VWICT = "Leadership beliefs about ICT outcomes", P_ICTUSE = "Principal ICT practices",
  P_ICTCOM = "Principal ICT communication", P_EXPTCH = "Teacher collaboration expectations",
  PV1CT = "Computational thinking achievement", PV2CT = "Computational thinking achievement",
  PV3CT = "Computational thinking achievement", PV4CT = "Computational thinking achievement",
  PV5CT = "Computational thinking achievement"
)
framework_role_map <- c(
  S_AGE = "Antecedent", S_NISB = "Antecedent", S_SPECEFF = "Process-related",
  S_GENEFF = "Process-related", S_USECOM = "Process-related", S_USEINF = "Process-related",
  S_ACCONT = "Process-related", S_SPECLASS = "Process-related", S_GENCLASS = "Process-related",
  S_CODLRN = "Process-related", S_ICTFUT = "Process-related", S_ICTLRN = "Process-related",
  S_ICTNEG = "Process-related", S_ICTPOS = "Process-related", S_GENACT = "Process-related",
  S_USESTD = "Process-related", PV1CIL = "Proximal achievement covariate",
  P_RATTCH = "Antecedent", P_EXPLRN = "Process-related", P_PRIORH = "Antecedent",
  P_PRIORS = "Antecedent", P_VWICT = "Process-related", P_ICTUSE = "Process-related",
  P_ICTCOM = "Process-related", P_EXPTCH = "Process-related",
  PV1CT = "Outcome plausible value", PV2CT = "Outcome plausible value",
  PV3CT = "Outcome plausible value", PV4CT = "Outcome plausible value",
  PV5CT = "Outcome plausible value"
)
variable_dictionary$construct <- unname(construct_map[variable_dictionary$variable])
variable_dictionary$framework_role <- unname(framework_role_map[variable_dictionary$variable])

kor_s <- read_sav(file.path(raw_dir, "BSGKORI2.sav"), n_max = 1)
kor_p <- read_sav(file.path(raw_dir, "BCGKORI2.sav"), n_max = 1)
variable_dictionary$label <- vapply(variable_dictionary$variable, function(v) {
  x <- if (v %in% names(kor_s)) kor_s[[v]] else kor_p[[v]]
  as.character(attr(x, "label"))
}, character(1))
write_csv_safe(variable_dictionary, "variable_dictionary.csv")

manifest <- data.frame(
  item = c(
    "analysis_date", "R_version", "seed", "official_data_source",
    "official_data_DOI", "primary_holdout", "cross_validation",
    "plausible_value_handling", "cil_pairing", "missing_data_sensitivity"
  ),
  value = c(
    as.character(Sys.Date()), R.version.string, seed_base,
    "https://www.iea.nl/data-tools/repository/icils",
    "https://doi.org/10.58150/ICILS_2018_data",
    "20% of schools held out; no student from a held-out school enters training",
    "Five folds assigned at the school level within training data",
    "Entire pipeline repeated for PV1CT-PV5CT; metrics averaged across plausible values",
    if (cil_mode == "matched") "CIL plausible value k paired with CT plausible value k (IEA convention)" else "PV1CIL used with every CT plausible value (fixed CIL pairing)",
    "Separate multiple_imputation.R: five predictor-only PMM imputations fitted using training rows only; all five CT PVs and two feature sets"
  )
)
write_csv_safe(manifest, "analysis_manifest.csv")

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Analysis complete. Results written to: ", normalizePath(output_dir))
