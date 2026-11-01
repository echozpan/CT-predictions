# Assess subgroup errors, calibration, and a screening scenario.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[[1]]
output_dir <- args[[2]]
source(file.path(script_dir, "ct_common.R"))
B <- 1000L

within_pred <- read.csv(file.path(output_dir, "within_country_predictions_all_plausible_values.csv"))
ext_pred <- read.csv(file.path(output_dir, "external_validation_predictions_all_plausible_values.csv"))
pred <- bind_rows(within_pred, ext_pred)

message("subgroup_audit: reading demographics...")
demo <- bind_rows(lapply(c("KOR", "DNK"), function(code) {
  d <- read_country(code)
  d[, c("ROW_UID", "S_SEX", "S_IMMIG", "S_NISB", "S_HOMLIT")]
}))
names(demo)[1] <- "row_uid"

student <- pred |>
  mutate(e = predicted - observed) |>
  group_by(context, target_code, feature_set, model, row_uid, school_uid) |>
  summarise(
    err = mean(e), sq = mean(e^2), abs_err = mean(abs(e)),
    observed = mean(observed), predicted = mean(predicted), .groups = "drop"
  ) |>
  left_join(demo, by = "row_uid")

base <- student |>
  group_by(context, target_code, feature_set, model) |>
  mutate(
    ses_tercile = cut(S_NISB, quantile(S_NISB, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                      labels = c("Lowest SES tercile", "Middle SES tercile", "Highest SES tercile"),
                      include.lowest = TRUE),
    sex = factor(S_SEX, levels = c(0, 1), labels = c("Boys", "Girls")),
    immigration = factor(ifelse(S_IMMIG == 0, "No immigrant background", "Immigrant background"),
                         levels = c("No immigrant background", "Immigrant background")),
    ct_quartile = cut(observed, quantile(observed, 0:4/4), labels = paste0("Observed CT Q", 1:4), include.lowest = TRUE)
  ) |>
  ungroup()

sub_stats <- function(d) {
  data.frame(
    n = nrow(d), mean_observed = mean(d$observed), mean_predicted = mean(d$predicted),
    mean_error = mean(d$err), RMSE = sqrt(mean(d$sq)), MAE = mean(d$abs_err)
  )
}
long <- bind_rows(
  base |> mutate(dimension = "Sex", subgroup = as.character(sex)),
  base |> mutate(dimension = "Socioeconomic tercile", subgroup = as.character(ses_tercile)),
  base |> mutate(dimension = "Immigration background", subgroup = as.character(immigration)),
  base |> mutate(dimension = "Observed CT quartile", subgroup = as.character(ct_quartile))
) |> filter(!is.na(subgroup))

by_sub <- long |>
  group_by(context, target_code, feature_set, model, dimension, subgroup) |>
  group_modify(~ sub_stats(.x)) |> ungroup()
write_csv_safe(by_sub, "fairness_subgroup_errors.csv")

contrast_def <- list(
  list(dim = "Sex", a = "Girls", b = "Boys"),
  list(dim = "Socioeconomic tercile", a = "Lowest SES tercile", b = "Highest SES tercile"),
  list(dim = "Immigration background", a = "Immigrant background", b = "No immigrant background")
)
boot_one <- function(d, dim, a, b, B) {
  d <- d[d$dimension == dim & d$subgroup %in% c(a, b), ]
  if (sum(d$subgroup == a) < 30 || sum(d$subgroup == b) < 30) return(NULL)
  schools <- unique(d$school_uid)
  idx_by_school <- split(seq_len(nrow(d)), d$school_uid)
  stat <- function(dd) {
    ea <- dd[dd$subgroup == a, ]; eb <- dd[dd$subgroup == b, ]
    c(bias_diff = mean(ea$err) - mean(eb$err),
      rmse_diff = sqrt(mean(ea$sq)) - sqrt(mean(eb$sq)))
  }
  est <- stat(d)
  set.seed(seed_base + 60000L)
  bs <- t(replicate(B, {
    pick <- sample(schools, length(schools), replace = TRUE)
    stat(d[unlist(idx_by_school[pick]), ])
  }))
  data.frame(
    dimension = dim, group_a = a, group_b = b,
    n_a = sum(d$subgroup == a), n_b = sum(d$subgroup == b),
    bias_diff = est[["bias_diff"]],
    bias_lo = quantile(bs[, 1], .025, na.rm = TRUE), bias_hi = quantile(bs[, 1], .975, na.rm = TRUE),
    rmse_diff = est[["rmse_diff"]],
    rmse_lo = quantile(bs[, 2], .025, na.rm = TRUE), rmse_hi = quantile(bs[, 2], .975, na.rm = TRUE),
    row.names = NULL
  )
}
keys <- long |> distinct(context, target_code, feature_set, model)
contrasts <- list()
for (i in seq_len(nrow(keys))) {
  kk <- keys[i, ]
  if (!(kk$model %in% c("Lasso", "Random forest"))) next
  dd <- long[long$context == kk$context & long$target_code == kk$target_code &
               long$feature_set == kk$feature_set & long$model == kk$model, ]
  for (cd in contrast_def) {
    r <- boot_one(dd, cd$dim, cd$a, cd$b, B)
    if (!is.null(r)) contrasts[[length(contrasts) + 1L]] <- cbind(kk, r)
  }
}
write_csv_safe(bind_rows(contrasts), "fairness_subgroup_contrasts_bootstrap.csv")

tr_cut <- list()
for (code in c("KOR", "DNK")) {
  d <- read_country(code)
  d <- d[complete.cases(d[, c(all_predictors, pv_names)]), ]
  flag <- make_group_holdout(d$SCHOOL_UID, prop = 0.20, seed = seed_base + match(code, c("KOR", "DNK")))
  tr_cut[[paste0("Within-country|", code)]] <- quantile(d$PV1CT[!flag], .25)
}
eu <- bind_rows(lapply(europe_train_codes, read_country))
eu <- eu[complete.cases(eu[, c(all_predictors, pv_names)]), ]
eu_cut <- unname(quantile(eu$PV1CT, .25))
screen <- pred |> filter(plausible_value == "PV1CT", model %in% c("Lasso", "Random forest")) |>
  left_join(demo, by = "row_uid") |>
  mutate(
    cut = ifelse(context == "Within-country",
                 unlist(tr_cut[paste0("Within-country|", target_code)]), eu_cut),
    truth = observed < cut, flag = predicted < cut,
    sex = factor(S_SEX, levels = c(0, 1), labels = c("Boys", "Girls"))
  )
screen_stats <- function(d) {
  tp <- sum(d$truth & d$flag); fn <- sum(d$truth & !d$flag)
  fp <- sum(!d$truth & d$flag); tn <- sum(!d$truth & !d$flag)
  data.frame(
    n = nrow(d), cut_score = d$cut[[1]], share_truly_below = mean(d$truth),
    share_flagged = mean(d$flag), sensitivity = tp / (tp + fn),
    false_negative_rate = fn / (tp + fn), precision = ifelse(tp + fp > 0, tp / (tp + fp), NA_real_),
    false_positive_rate = fp / (fp + tn), accuracy = (tp + tn) / nrow(d)
  )
}
screen_all <- screen |> group_by(context, target_code, feature_set, model) |>
  group_modify(~ screen_stats(.x)) |> ungroup() |> mutate(subgroup = "All students")
screen_sex <- screen |> filter(!is.na(sex)) |> group_by(context, target_code, feature_set, model, subgroup = as.character(sex)) |>
  group_modify(~ screen_stats(.x)) |> ungroup()
write_csv_safe(bind_rows(screen_all, screen_sex), "fairness_screening_scenario.csv")

sens_boot <- list()
sk <- screen |> filter(!is.na(sex)) |> distinct(context, target_code, feature_set, model)
for (i in seq_len(nrow(sk))) {
  kk <- sk[i, ]
  dd <- screen[screen$context == kk$context & screen$target_code == kk$target_code &
                 screen$feature_set == kk$feature_set & screen$model == kk$model & !is.na(screen$sex), ]
  sens <- function(x, g) { z <- x[x$sex == g & x$truth, ]; mean(z$flag) }
  est <- sens(dd, "Girls") - sens(dd, "Boys")
  idx <- split(seq_len(nrow(dd)), dd$school_uid)
  schools <- names(idx)
  set.seed(seed_base + 61000L + i)
  bs <- replicate(B, {
    pick <- sample(schools, length(schools), replace = TRUE)
    x <- dd[unlist(idx[pick]), ]
    sens(x, "Girls") - sens(x, "Boys")
  })
  sens_boot[[i]] <- cbind(kk, data.frame(
    n_truly_below_girls = sum(dd$sex == "Girls" & dd$truth), n_truly_below_boys = sum(dd$sex == "Boys" & dd$truth),
    sensitivity_girls = sens(dd, "Girls"), sensitivity_boys = sens(dd, "Boys"),
    sensitivity_diff_girls_minus_boys = est,
    lo = unname(quantile(bs, .025, na.rm = TRUE)), hi = unname(quantile(bs, .975, na.rm = TRUE)), row.names = NULL))
}
write_csv_safe(bind_rows(sens_boot), "fairness_screening_sex_sensitivity_bootstrap.csv")

sex_desc <- bind_rows(lapply(c("KOR", "DNK"), function(code) {
  d <- read_country(code)
  d <- d[complete.cases(d[, c(all_predictors, pv_names)]) & d$S_SEX %in% c(0, 1), ]
  d$sex <- ifelse(d$S_SEX == 1, "Girls", "Boys")
  d |> group_by(sex) |> summarise(n = n(), mean_PV1CIL = mean(PV1CIL), mean_PV1CT = mean(PV1CT),
                                  .groups = "drop") |> mutate(target_code = code, .before = 1)
}))
write_csv_safe(sex_desc, "fairness_sex_cil_ct_descriptives.csv")
message("subgroup_audit complete.")
