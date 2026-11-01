# Compare alternative country compositions of the training pool.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "ct_common.R"))

message("training_pool_composition: reading countries...")
country_data <- lapply(names(country_map), read_country)
names(country_data) <- names(country_map)
cc <- lapply(country_data, function(d) d[complete.cases(d[, c(all_predictors, pv_names)]), ])

holdouts <- list()
for (code in c("KOR", "DNK")) {
  d <- cc[[code]]
  flag <- make_group_holdout(d$SCHOOL_UID, prop = 0.20, seed = seed_base + match(code, c("KOR", "DNK")))
  holdouts[[code]] <- d[flag, ]
}

pools <- list(
  "Europe-5 (reference pool)" = europe_train_codes,
  "Core-3 (Finland, France, Germany)" = c("FIN", "FRA", "DEU"),
  "Europe-5 without Finland" = setdiff(europe_train_codes, "FIN"),
  "Europe-5 without France" = setdiff(europe_train_codes, "FRA"),
  "Europe-5 without Germany" = setdiff(europe_train_codes, "DEU"),
  "Europe-5 without Luxembourg" = setdiff(europe_train_codes, "LUX"),
  "Europe-5 without Portugal" = setdiff(europe_train_codes, "PRT"),
  "Europe-5 plus United States" = c(europe_train_codes, "USA"),
  "United States only" = "USA",
  "Denmark only (target Korea)" = "DNK",
  "Korea only (target Denmark)" = "KOR"
)
spec <- bind_rows(
  expand.grid(model = "Lasso", pv = pv_names, stringsAsFactors = FALSE),
  data.frame(model = "Random forest", pv = "PV1CT", stringsAsFactors = FALSE)
)
feature_sets <- list(
  "All predictors" = all_predictors,
  "Excluding PV1CIL" = setdiff(all_predictors, "PV1CIL")
)

rows <- list()
k <- 0L
for (pool_name in names(pools)) {
  members <- pools[[pool_name]]
  ds <- bind_rows(cc[members])
  targets <- c("KOR", "DNK")
  targets <- targets[!(targets %in% members)]
  for (fs in names(feature_sets)) {
    pred_vars <- feature_sets[[fs]]
    for (i in seq_len(nrow(spec))) {
      k <- k + 1L
      model_name <- spec$model[[i]]
      pv <- spec$pv[[i]]
      message("  pool ", pool_name, " | ", fs, " | ", pv, " | ", model_name)
      fit <- fit_model(model_name, pred_matrix(ds, pred_vars, pv), ds[[pv]],
                       ds$SCHOOL_UID, seed_base + 40000L + k)
      for (tg in targets) {
        dt <- holdouts[[tg]]
        prd <- predict_model(fit, pred_matrix(dt, pred_vars, pv))
        rows[[length(rows) + 1L]] <- data.frame(
          pool = pool_name, pool_members = paste(members, collapse = "+"),
          target = tg, feature_set = fs, model = model_name,
          plausible_value = pv, train_students = nrow(ds),
          train_schools = n_distinct(ds$SCHOOL_UID),
          train_mean_CT = mean(ds[[pv]]), train_sd_CT = sd(ds[[pv]]),
          test_students = nrow(dt), t(metric_values(dt[[pv]], prd)),
          row.names = NULL
        )
      }
    }
  }
}
res <- bind_rows(rows)
write_csv_safe(res, "training_pool_composition_each_PV.csv")
pooled <- res |>
  group_by(pool, pool_members, target, feature_set, model) |>
  summarise(
    plausible_values = n(), train_students = first(train_students),
    train_schools = first(train_schools), test_students = first(test_students),
    RMSE = mean(RMSE), MAE = mean(MAE), R2 = mean(R2),
    mean_error = mean(mean_error), RMSE_over_target_SD = mean(RMSE_over_target_SD),
    .groups = "drop"
  )
write_csv_safe(pooled, "training_pool_composition_pooled.csv")
message("training_pool_composition complete: ", nrow(res), " rows")
