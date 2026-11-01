# Test associations between transfer penalties and individualism distance.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[[1]]
output_dir <- args[[2]]
source(file.path(script_dir, "ct_common.R"))

idv <- c(DEU = 67, DNK = 74, FIN = 63, FRA = 71, KOR = 18, LUX = 60, PRT = 27, USA = 91)
codes <- names(idv)

pen <- read.csv(file.path(output_dir, "transfer_penalties.csv"))
pen <- pen[pen$source %in% codes, ]

cc <- lapply(codes, function(code) {
  d <- read_country(code); d[complete.cases(d[, c(all_predictors, pv_names)]), ]
})
names(cc) <- codes
smd <- function(a, b, vars) {
  mean(sapply(vars, function(v) {
    pooled <- sqrt((var(a[[v]]) + var(b[[v]])) / 2)
    abs(mean(b[[v]]) - mean(a[[v]])) / pooled
  }))
}
pair_shift <- expand.grid(source = codes, target = codes, stringsAsFactors = FALSE)
pair_shift <- pair_shift[pair_shift$source != pair_shift$target, ]
pair_shift$smd_all <- mapply(function(s, t) smd(cc[[s]], cc[[t]], all_predictors), pair_shift$source, pair_shift$target)
pair_shift$smd_noCIL <- mapply(function(s, t) smd(cc[[s]], cc[[t]], setdiff(all_predictors, "PV1CIL")), pair_shift$source, pair_shift$target)
write_csv_safe(pair_shift, "pairwise_predictor_shift.csv")

dyads <- pen |>
  filter(model == "Lasso") |>
  left_join(pair_shift, by = c("source", "target")) |>
  mutate(
    idv_distance = abs(idv[source] - idv[target]),
    idv_signed = idv[target] - idv[source],
    ct_mean_distance = abs(mean_shift_CT),
    ct_logsd_distance = abs(log(sd_ratio_CT)),
    shift = ifelse(feature_set == "All predictors", smd_all, smd_noCIL)
  )
write_csv_safe(dyads, "cultural_distance_dyads.csv")

perms <- function(v) {
  if (length(v) <= 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) for (p in perms(v[-i])) out[[length(out) + 1L]] <- c(v[i], p)
  out
}
all_perm <- perms(seq_along(codes))

res <- list()
for (fs in unique(dyads$feature_set)) {
  d <- dyads[dyads$feature_set == fs, ]
  d$si <- match(d$source, codes); d$ti <- match(d$target, codes)
  stat_sp <- function(idv_vec) {
    dist <- abs(idv_vec[d$si] - idv_vec[d$ti])
    suppressWarnings(cor(dist, d$transfer_penalty_RMSE, method = "spearman"))
  }
  stat_beta <- function(idv_vec) {
    dist <- abs(idv_vec[d$si] - idv_vec[d$ti])
    coef(lm(d$transfer_penalty_RMSE ~ dist + d$shift + d$ct_mean_distance + d$ct_logsd_distance))[2]
  }
  obs_sp <- stat_sp(idv); obs_beta <- stat_beta(idv)
  null_sp <- sapply(all_perm, function(p) stat_sp(idv[p]))
  null_beta <- sapply(all_perm, function(p) stat_beta(idv[p]))
  fit <- lm(transfer_penalty_RMSE ~ idv_distance + shift + ct_mean_distance + ct_logsd_distance, data = d)
  res[[length(res) + 1L]] <- data.frame(
    feature_set = fs, dyads = nrow(d),
    spearman_idv_distance_penalty = obs_sp,
    permutation_p_spearman_two_sided = mean(abs(null_sp) >= abs(obs_sp) - 1e-12, na.rm = TRUE),
    adj_beta_idv_per_10_points = 10 * unname(obs_beta),
    permutation_p_adj_beta_two_sided = mean(abs(null_beta) >= abs(obs_beta) - 1e-12, na.rm = TRUE),
    spearman_predictor_shift_penalty = cor(d$shift, d$transfer_penalty_RMSE, method = "spearman"),
    spearman_ct_mean_distance_penalty = cor(d$ct_mean_distance, d$transfer_penalty_RMSE, method = "spearman"),
    spearman_ct_sd_ratio_penalty = cor(d$ct_logsd_distance, d$transfer_penalty_RMSE, method = "spearman"),
    r2_multivariable = summary(fit)$r.squared
  )
}
write_csv_safe(bind_rows(res), "cultural_distance_mantel_tests.csv")

tgt <- dyads |> group_by(feature_set, target) |>
  summarise(mean_penalty = mean(transfer_penalty_RMSE), mean_R2_loss = mean(transfer_penalty_R2), .groups = "drop") |>
  mutate(idv = idv[target])
write_csv_safe(tgt, "target_level_penalty_vs_idv.csv")
tl <- bind_rows(lapply(unique(tgt$feature_set), function(fs) {
  t <- tgt[tgt$feature_set == fs, ]
  ct <- cor.test(t$idv, t$mean_penalty, method = "spearman", exact = TRUE)
  data.frame(feature_set = fs, n_targets = nrow(t), spearman_rho = unname(ct$estimate), exact_p = ct$p.value)
}))
write_csv_safe(tl, "target_level_penalty_vs_idv_test.csv")
message("cultural_distance complete.")
