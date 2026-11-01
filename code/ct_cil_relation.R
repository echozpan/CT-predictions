# Summarize CIL–CT associations and offsets from the European training pool.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[[1]]
output_dir <- args[[2]]
source(file.path(script_dir, "ct_common.R"))

country_data <- lapply(names(country_map), read_country)
names(country_data) <- names(country_map)
cc <- lapply(country_data, function(d) d[complete.cases(d[, c(all_predictors, pv_names)]), ])
eur <- bind_rows(cc[europe_train_codes])

rows <- list()
for (code in names(cc)) {
  d <- cc[[code]]
  per_pv <- lapply(pv_names, function(pv) {
    cil <- d[[cil_col_for(pv)]]; ct <- d[[pv]]
    fit_eur <- lm(eur[[pv]] ~ eur[[cil_col_for(pv)]])
    pred <- coef(fit_eur)[1] + coef(fit_eur)[2] * cil
    data.frame(r = cor(cil, ct), slope = unname(coef(lm(ct ~ cil))[2]),
               mean_CIL = mean(cil), mean_CT = mean(ct), sd_CIL = sd(cil), sd_CT = sd(ct),
               ct_resid_vs_europe = mean(ct - pred))
  })
  m <- colMeans(bind_rows(per_pv))
  rows[[code]] <- data.frame(country_code = code, country = unname(country_map[[code]]), n = nrow(d), t(m), row.names = NULL)
}
res <- bind_rows(rows)
res$R2_CIL_only_implied <- res$r^2
write_csv_safe(res, "ct_cil_relation_by_country.csv")
print(res, digits = 3)
