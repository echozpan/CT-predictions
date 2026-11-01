# Compare cross-country transfer with school-grouped local benchmarks.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "ct_common.R"))

message("transfer_matrix: reading countries...")
country_data <- lapply(names(country_map), read_country)
names(country_data) <- names(country_map)
cc <- lapply(country_data, function(d) d[complete.cases(d[, c(all_predictors, pv_names)]), ])

sources <- c(as.list(cc), list(EUR5 = bind_rows(cc[europe_train_codes])))
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
for (src in names(sources)) {
  ds <- sources[[src]]
  targets <- if (src == "EUR5") c("KOR", "DNK", "USA") else names(cc)
  for (fs in names(feature_sets)) {
    pred_vars <- feature_sets[[fs]]
    for (i in seq_len(nrow(spec))) {
      model_name <- spec$model[[i]]
      pv <- spec$pv[[i]]
      k <- k + 1L
      message("  source ", src, " | ", fs, " | ", pv, " | ", model_name)
      x <- pred_matrix(ds, pred_vars, pv)
      y <- ds[[pv]]
      g <- ds$SCHOOL_UID
      fit <- fit_model(model_name, x, y, g, seed_base + 30000L + k)
      for (tg in targets) {
        dt <- cc[[tg]]
        if (tg == src) {

          fold <- make_group_foldid(g, k = 5L, seed = seed_base + 31000L + k)
          oof <- rep(NA_real_, nrow(ds))
          for (f in sort(unique(fold))) {
            tr <- fold != f
            fo <- fit_model(model_name, x[tr, , drop = FALSE], y[tr], g[tr],
                            seed_base + 32000L + 10L * k + f)
            oof[!tr] <- predict_model(fo, x[!tr, , drop = FALSE])
          }
          obs <- y
          prd <- oof
        } else {
          obs <- dt[[pv]]
          prd <- predict_model(fit, pred_matrix(dt, pred_vars, pv))
        }
        m <- metric_values(obs, prd)
        rows[[length(rows) + 1L]] <- data.frame(
          source = src, target = tg, is_within = (src == tg), feature_set = fs,
          model = model_name, plausible_value = pv,
          train_students = nrow(ds), target_students = length(obs),
          t(m), correlation = cor(obs, prd), source_mean_CT = mean(y),
          source_sd_CT = sd(y), row.names = NULL
        )
      }
    }
  }
}
res <- bind_rows(rows)
write_csv_safe(res, "transfer_matrix_all_pairs.csv")

pooled <- res |>
  group_by(source, target, is_within, feature_set, model) |>
  summarise(
    plausible_values = n(), train_students = first(train_students),
    target_students = first(target_students),
    RMSE = mean(RMSE), MAE = mean(MAE), R2 = mean(R2),
    mean_error = mean(mean_error), correlation = mean(correlation),
    target_mean_CT = mean(target_mean), target_sd_CT = mean(target_sd),
    source_mean_CT = mean(source_mean_CT), source_sd_CT = mean(source_sd_CT),
    .groups = "drop"
  )
write_csv_safe(pooled, "transfer_matrix_pooled.csv")

within <- pooled |> filter(is_within) |>
  select(target, feature_set, model, RMSE_within = RMSE, R2_within = R2)
pen <- pooled |> filter(!is_within) |>
  left_join(within, by = c("target", "feature_set", "model")) |>
  mutate(
    transfer_penalty_RMSE = RMSE - RMSE_within,
    transfer_penalty_R2 = R2 - R2_within,
    relative_penalty_percent = 100 * (RMSE / RMSE_within - 1),
    mean_shift_CT = target_mean_CT - source_mean_CT,
    sd_ratio_CT = target_sd_CT / source_sd_CT
  )
write_csv_safe(pen, "transfer_penalties.csv")
message("transfer_matrix complete: ", nrow(res), " metric rows")
