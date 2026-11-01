# Evaluate fixed predictions using student sampling weights.
suppressPackageStartupMessages({library(dplyr);library(haven)})
a <- commandArgs(TRUE);raw <- a[1];out <- a[2]
p <- bind_rows(read.csv(file.path(out,'within_country_predictions_all_plausible_values.csv')),
               read.csv(file.path(out,'external_validation_predictions_all_plausible_values.csv')))
ws <- bind_rows(lapply(c('KOR','DNK'),function(code){
  s <- read_sav(file.path(raw,paste0('BSG',code,'I2.sav')),
    col_select=all_of(c('IDSTUD','TOTWGTS')))
  data.frame(row_uid=paste(code,s$IDSTUD,sep='_'),weight=as.numeric(s$TOTWGTS))
}))
stopifnot(!anyDuplicated(ws$row_uid))
p <- left_join(p,ws,by='row_uid')
stopifnot(all(is.finite(p$weight)),all(p$weight>0))
r <- p |> group_by(context,target_country,target_code,feature_set,plausible_value,model) |>
  summarise(n=n(),kish_effective_n=sum(weight)^2/sum(weight^2),
    MSE=weighted.mean((predicted-observed)^2,weight),RMSE=sqrt(MSE),
    MAE=weighted.mean(abs(predicted-observed),weight),
    R2=1-sum(weight*(predicted-observed)^2)/sum(weight*(observed-weighted.mean(observed,weight))^2),
    mean_error=weighted.mean(predicted-observed,weight),.groups='drop')
pooled <- r |> group_by(context,target_country,target_code,feature_set,model) |>
  summarise(plausible_values=n(),n=mean(n),kish_effective_n=mean(kish_effective_n),
    RMSE_mean=mean(RMSE),RMSE_SD_across_PVs=sd(RMSE),MAE_mean=mean(MAE),
    R2_mean=mean(R2),mean_error_mean=mean(mean_error),.groups='drop')
write.csv(r,file.path(out,'weighted_evaluation_metrics_each_PV.csv'),row.names=FALSE)
write.csv(pooled,file.path(out,'weighted_evaluation_metrics_pooled.csv'),row.names=FALSE)
