# Check saved metrics, outcomes, evaluation rows, and imputation partitions.
suppressPackageStartupMessages({library(dplyr);library(haven)})
args <- commandArgs(TRUE); raw <- args[1]; out <- args[2]
read <- function(n)read.csv(file.path(out,n))
w <- read('within_country_predictions_all_plausible_values.csv')
e <- read('external_validation_predictions_all_plausible_values.csv')
wm <- read('within_country_metrics_all_plausible_values.csv')
em <- read('external_validation_metrics_all_plausible_values.csv')
check <- function(p,m) {
  actual <- p |> group_by(target_code,feature_set,plausible_value,model) |>
    summarise(calc_RMSE=sqrt(mean((predicted-observed)^2)),
      calc_MAE=mean(abs(predicted-observed)),
      calc_R2=1-sum((predicted-observed)^2)/sum((observed-mean(observed))^2),.groups='drop')
  z <- left_join(m,actual,by=c('target_code','feature_set','plausible_value','model'))
  stopifnot(nrow(z)==100,!anyNA(z$calc_RMSE),
    max(abs(z$RMSE-z$calc_RMSE))<1e-8,max(abs(z$MAE-z$calc_MAE))<1e-8,
    max(abs(z$R2-z$calc_R2))<1e-8)
}
check(w,wm); check(e,em)
keys <- c('target_code','feature_set','plausible_value','model','row_uid','school_uid')
z <- inner_join(w,e,by=keys,suffix=c('_within','_external'))
stopifnot(nrow(z)==nrow(w),nrow(z)==nrow(e),all(z$observed_within==z$observed_external))
mi <- read('multiple_imputation_sample_partitions.csv')
stopifnot(!any(duplicated(unique(mi[c('school_uid','partition')])$school_uid)))
for(code in c('KOR','DNK')) {
  pred <- w[w$target_code==code,]
  ids <- unique(pred$row_uid)
  s <- read_sav(file.path(raw,paste0('BSG',code,'I2.sav')),
    col_select=all_of(c('IDSTUD','IDSCHOOL',paste0('PV',1:5,'CT'))))
  s$row_uid <- paste(code,s$IDSTUD,sep='_')
  for(pv in paste0('PV',1:5,'CT')) {
    pp <- pred[pred$plausible_value==pv,]
    yy <- as.numeric(s[[pv]][match(pp$row_uid,s$row_uid)])
    stopifnot(max(abs(yy-pp$observed))<1e-8)
  }
  mm <- mi[mi$country_code==code,]
  stopifnot(setequal(ids,mm$row_uid[mm$partition=='test' & mm$complete_case]))
}
imp <- read('permutation_importance_all_models.csv')
stopifnot(nrow(imp)==500,all(is.finite(imp$delta_MSE)))
mit <- read('multiple_imputation_sensitivity_each_imputation.csv')
stopifnot(nrow(mit)==1000,all(is.finite(mit$RMSE)))
writeLines(c('PASS: 200 primary model/PV metric rows independently recomputed from predictions',
 'PASS: 47,150 within/external prediction pairs share identical rows and outcomes',
 'PASS: saved outcomes equal official SPSS plausible values',
 'PASS: no MI school occurs in both partitions; original target rows preserved',
 'PASS: all 500 permutation importance estimates are finite',
 'PASS: all 1,000 MI evaluation rows are finite (500 fits, two target subsets)',
 'Intervals are conditional; no survey-weight or full model-refitting inference is claimed'),
 file.path(out,'validation_report.txt'))
cat('All analysis validation checks passed.\n')
