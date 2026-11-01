# Shared data preparation, school splits, model fitting, and evaluation functions.
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

options(stringsAsFactors = FALSE, scipen = 999)

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
models <- c("Lasso", "Ridge", "PCR", "Random forest", "KNN")
