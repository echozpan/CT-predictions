# Compare theoretical predictor blocks alone and with each block removed.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "ct_common.R"))

blocks <- list(
  "Student antecedents" = c("S_AGE", "S_NISB"),
  "Efficacy and attitudes" = c("S_SPECEFF", "S_GENEFF", "S_ICTPOS", "S_ICTNEG", "S_ICTFUT"),
  "ICT use" = c("S_USECOM", "S_USEINF", "S_ACCONT", "S_GENACT", "S_USESTD"),
  "Learning opportunities" = c("S_SPECLASS", "S_GENCLASS", "S_CODLRN", "S_ICTLRN"),
  "School context" = school_predictors
)
stopifnot(setequal(c(unlist(blocks), "PV1CIL"), all_predictors))
no_cil <- setdiff(all_predictors, "PV1CIL")

sets <- list(
  "All predictors" = all_predictors,
  "Excluding PV1CIL" = no_cil,
  "CIL only" = "PV1CIL",
  "Student-level only (with CIL)" = student_predictors,
  "Student-level only (no CIL)" = setdiff(student_predictors, "PV1CIL"),
  "School-level only" = school_predictors
)
for (b in names(blocks)) {
  sets[[paste0("No-CIL minus ", b)]] <- setdiff(no_cil, blocks[[b]])
  sets[[paste0("No-CIL ", b, " only")]] <- blocks[[b]]
}

message("feature_block_ablation: reading countries...")
country_data <- lapply(names(country_map), read_country)
names(country_data) <- names(country_map)
cc <- lapply(country_data, function(d) d[complete.cases(d[, c(all_predictors, pv_names)]), ])

holdout <- list()
for (code in c("KOR", "DNK")) {
  d <- cc[[code]]
  flag <- make_group_holdout(d$SCHOOL_UID, prop = 0.20, seed = seed_base + match(code, c("KOR", "DNK")))
  holdout[[code]] <- list(train = d[!flag, ], test = d[flag, ])
}
europe <- bind_rows(cc[europe_train_codes])

mk <- function(d, vars, pv) {
  m <- pred_matrix(d, vars, pv)
  if (ncol(m) == 1L) m <- cbind(m, m)
  m
}

spec <- bind_rows(
  expand.grid(model = "Lasso", pv = pv_names, stringsAsFactors = FALSE),
  data.frame(model = "Random forest", pv = "PV1CT", stringsAsFactors = FALSE)
)

rows <- list()
k <- 0L
for (set_name in names(sets)) {
  pv_set <- sets[[set_name]]
  for (i in seq_len(nrow(spec))) {
    model_name <- spec$model[[i]]
    pv <- spec$pv[[i]]

    for (code in c("KOR", "DNK")) {
      k <- k + 1L
      message("  ", set_name, " | within ", code, " | ", pv, " | ", model_name)
      tr <- holdout[[code]]$train
      te <- holdout[[code]]$test
      fit <- fit_model(model_name, mk(tr, pv_set, pv), tr[[pv]], tr$SCHOOL_UID, seed_base + 50000L + k)
      prd <- predict_model(fit, mk(te, pv_set, pv))
      rows[[length(rows) + 1L]] <- data.frame(
        context = "Within-country", target = code, feature_set = set_name,
        n_predictors = length(pv_set), model = model_name, plausible_value = pv,
        t(metric_values(te[[pv]], prd)), row.names = NULL
      )
    }

    k <- k + 1L
    message("  ", set_name, " | Europe-trained | ", pv, " | ", model_name)
    fit <- fit_model(model_name, mk(europe, pv_set, pv), europe[[pv]], europe$SCHOOL_UID, seed_base + 50000L + k)
    for (code in c("KOR", "DNK")) {
      te <- holdout[[code]]$test
      prd <- predict_model(fit, mk(te, pv_set, pv))
      rows[[length(rows) + 1L]] <- data.frame(
        context = "Europe-trained external", target = code, feature_set = set_name,
        n_predictors = length(pv_set), model = model_name, plausible_value = pv,
        t(metric_values(te[[pv]], prd)), row.names = NULL
      )
    }
  }
}
res <- bind_rows(rows)
write_csv_safe(res, "feature_block_ablation_each_PV.csv")
pooled <- res |>
  group_by(context, target, feature_set, n_predictors, model) |>
  summarise(plausible_values = n(), RMSE = mean(RMSE), MAE = mean(MAE), R2 = mean(R2),
            RMSE_over_target_SD = mean(RMSE_over_target_SD), .groups = "drop")
write_csv_safe(pooled, "feature_block_ablation_pooled.csv")
message("feature_block_ablation complete: ", nrow(res), " rows")
