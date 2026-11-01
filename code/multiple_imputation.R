# Evaluate prediction after training-only imputation of missing predictors.
suppressPackageStartupMessages({
  library(dplyr); library(haven); library(mice); library(glmnet)
  library(RANN); library(ranger)
})
args <- commandArgs(TRUE)
stopifnot(length(args) >= 2)
raw_dir <- args[1]; out <- args[2]
seed_base <- 20260812L
cil_mode <- Sys.getenv('CIL_MODE','matched')
stopifnot(cil_mode %in% c('matched','fixed'))
cil_col_for <- function(pv) paste0('PV',sub('^PV([0-9])CT$','\\1',pv),'CIL')

xpv <- function(dd,vars,pv,rows) {
  x <- dd[rows,vars]
  if(cil_mode=='matched' && 'PV1CIL' %in% vars) x$PV1CIL <- dd[[cil_col_for(pv)]][rows]
  x
}
student <- c('S_AGE','S_SPECEFF','S_GENEFF','S_USECOM','S_USEINF',
 'S_ACCONT','S_SPECLASS','S_GENCLASS','S_CODLRN','S_ICTFUT','S_ICTLRN',
 'S_ICTNEG','S_ICTPOS','S_GENACT','S_USESTD','S_NISB','PV1CIL')
school <- c('P_RATTCH','P_EXPLRN','P_PRIORH','P_PRIORS','P_VWICT',
 'P_ICTUSE','P_ICTCOM','P_EXPTCH')
aux <- c('S_SEX','S_HISEI','S_HOMLIT','S_IMMIG')
pvs <- paste0('PV',1:5,'CT'); models <- c('Lasso','Ridge','PCR','Random forest','KNN')
only_pcr <- length(args) >= 3 && args[3] == 'PCR'
if(only_pcr) models <- 'PCR'
hp <- read.csv(file.path(out,'selected_hyperparameters.csv'))
wp <- read.csv(file.path(out,'within_country_predictions_all_plausible_values.csv'))
ccm <- read.csv(file.path(out,'within_country_metrics_all_plausible_values.csv'))
rows <- preds <- split_rows <- logs <- list(); counter <- 0L

pmm <- function(dat, test, m, seed) {
  stopifnot(!any(pvs %in% names(dat)), sum(!test)>10)
  init <- dat
  for (v in names(init)) {
    donor <- dat[[v]][!test & is.finite(dat[[v]])]
    if (!length(donor)) stop('No observed training donor for ',v)
    init[[v]][is.na(init[[v]])] <- median(donor)
  }
  fit <- mice(dat, m=m, maxit=10, method='pmm', ignore=test,
    data.init=init, remove.collinear=FALSE, remove.constant=FALSE,
    printFlag=FALSE, seed=seed)
  if (!is.null(fit$loggedEvents)) {
    ev <- fit$loggedEvents; ev$seed <- seed
    logs[[length(logs)+1L]] <<- ev
  }
  lapply(seq_len(m), function(i) {
    z <- complete(fit,i); stopifnot(!anyNA(z)); z
  })
}
scale_pair <- function(a,b) {
  mu <- colMeans(a); sig <- apply(a,2,sd); sig[!is.finite(sig)|sig==0] <- 1
  list(a=sweep(sweep(a,2,mu,'-'),2,sig,'/'),
       b=sweep(sweep(b,2,mu,'-'),2,sig,'/'))
}
predict_fixed <- function(model,a,y,b,h,seed) {
  a <- as.matrix(a); b <- as.matrix(b)
  if(model %in% c('Lasso','Ridge')) {
    f <- glmnet(a,y,alpha=ifelse(model=='Lasso',1,0),lambda=h,standardize=TRUE)
    return(as.numeric(predict(f,newx=b,s=h)))
  }
  if(model=='Random forest') {
    f <- ranger(x=a,y=y,mtry=as.integer(h),num.trees=500,
      min.node.size=5,num.threads=1,seed=seed)
    return(as.numeric(predict(f,data=b)$predictions))
  }
  sc <- scale_pair(a,b)
  if(model=='PCR') {
    nc <- as.integer(h); pc <- prcomp(sc$a,center=FALSE,scale.=FALSE,rank.=nc)
    co <- qr.solve(cbind(1,pc$x[,seq_len(nc),drop=FALSE]),y)
    return(as.numeric(cbind(1,sc$b%*%pc$rotation[,seq_len(nc),drop=FALSE])%*%co))
  }
  ix <- RANN::nn2(sc$a,sc$b,k=as.integer(h))$nn.idx
  rowMeans(matrix(y[ix],nrow=nrow(ix)))
}
metric <- function(y,p) {
  e <- p-y; mse <- mean(e^2)
  c(n=length(y),MSE=mse,RMSE=sqrt(mse),MAE=mean(abs(e)),
    R2=1-sum(e^2)/sum((y-mean(y))^2),mean_error=mean(e),
    target_mean=mean(y),target_sd=sd(y),RMSE_in_ICILS_SD=sqrt(mse)/100,
    RMSE_over_target_SD=sqrt(mse)/sd(y))
}
for(code in c('KOR','DNK')) {
  country <- ifelse(code=='KOR','Korea','Denmark')
  s <- read_sav(file.path(raw_dir,paste0('BSG',code,'I2.sav')),
    col_select=all_of(c('IDSCHOOL','IDSTUD',student,aux,pvs,paste0('PV',2:5,'CIL'))))
  q <- read_sav(file.path(raw_dir,paste0('BCG',code,'I2.sav')),
    col_select=all_of(c('IDSCHOOL',school)))
  stopifnot(!anyDuplicated(q$IDSCHOOL))
  d <- left_join(s,q,by='IDSCHOOL'); d[] <- lapply(d,function(x)as.numeric(zap_labels(x)))
  d <- d[complete.cases(d[,pvs]),]
  d$school_uid <- paste(code,d$IDSCHOOL,sep='_')
  d$row_uid <- paste(code,d$IDSTUD,sep='_')
  cc <- complete.cases(d[,c(student,school,pvs)])
  primary_test <- unique(wp$school_uid[wp$target_code==code])
  cc_schools <- unique(d$school_uid[cc])
  extra <- setdiff(unique(d$school_uid),cc_schools)
  set.seed(seed_base+match(code,c('KOR','DNK')))
  extra_test <- if(length(extra))sample(extra,round(.2*length(extra))) else character()
  test <- d$school_uid %in% c(primary_test,extra_test)
  stopifnot(!length(intersect(d$school_uid[test],d$school_uid[!test])))
  stopifnot(setequal(d$row_uid[test & cc],unique(wp$row_uid[wp$target_code==code])))
  split_rows[[code]] <- data.frame(country_code=code,row_uid=d$row_uid,
    school_uid=d$school_uid,partition=ifelse(test,'test','train'),complete_case=cc)
  sch <- d[!duplicated(d$school_uid),c('school_uid',school)]
  sch_test <- sch$school_uid %in% unique(d$school_uid[test])
  school_sets <- pmm(sch[,school],sch_test,5,seed_base+31000+match(code,c('KOR','DNK')))
  for(fs in c('All predictors','Excluding PV1CIL')) {
    st <- if(fs=='All predictors')student else setdiff(student,'PV1CIL')
    vars <- c(st,school)
    for(i in 1:5) {
      dd <- d; dd[,school] <- school_sets[[i]][match(d$school_uid,sch$school_uid),]

      si <- pmm(dd[,c(st,school,aux)],test,1,
        seed_base+32000+100*match(code,c('KOR','DNK'))+10*match(fs,c('All predictors','Excluding PV1CIL'))+i)[[1]]
      dd[,names(si)] <- si
      for(pv in pvs) for(model in models) {
        counter <- counter+1L
        h <- hp$hyperparameter[hp$context=='Within-country' & hp$target_code==code &
          hp$feature_set==fs & hp$plausible_value==pv & hp$model==model]
        stopifnot(length(h)==1,is.finite(h))
        message(code,' ',fs,' imputation ',i,' ',pv,' ',model)
        pr <- predict_fixed(model,xpv(dd,vars,pv,!test),dd[[pv]][!test],xpv(dd,vars,pv,test),h,
          seed_base+40000+counter)
        for(subset in c('All eligible held-out students','Complete-case held-out subset')) {
          keep <- if(subset=='All eligible held-out students')rep(TRUE,sum(test)) else cc[test]
          rows[[length(rows)+1L]] <- cbind(data.frame(country=country,country_code=code,
            feature_set=fs,plausible_value=pv,imputation=i,model=model,
            evaluation_subset=subset,train_n=sum(!test),test_n=sum(keep),hyperparameter=h),
            as.data.frame(as.list(metric(dd[[pv]][test][keep],pr[keep]))))
        }
        preds[[counter]] <- data.frame(country_code=code,feature_set=fs,plausible_value=pv,
          imputation=i,model=model,row_uid=d$row_uid[test],school_uid=d$school_uid[test],
          complete_case=cc[test],observed=d[[pv]][test],predicted=pr)
      }
    }
  }
}
each <- bind_rows(rows)
if(only_pcr) {
 old <- read.csv(file.path(out,'multiple_imputation_sensitivity_each_imputation.csv'))
 each <- bind_rows(old[old$model != 'PCR',],each)
 prediction_file <- file.path(out,'multiple_imputation_predictions.csv')
 if (!file.exists(prediction_file)) prediction_file <- gzfile(paste0(prediction_file,'.gz'))
 oldp <- read.csv(prediction_file)
 preds <- c(list(oldp[oldp$model != 'PCR',]),preds)
}
pooled <- each |> group_by(country,country_code,feature_set,model,evaluation_subset) |>
  summarise(imputations=n_distinct(imputation),plausible_values=n_distinct(plausible_value),
    train_n_mean=mean(train_n),test_n_mean=mean(test_n),RMSE_mean=mean(RMSE),
    RMSE_SD=sd(RMSE),MAE_mean=mean(MAE),MAE_SD=sd(MAE),R2_mean=mean(R2),R2_SD=sd(R2),
    RMSE_in_ICILS_SD_mean=mean(RMSE_in_ICILS_SD),.groups='drop')
write.csv(each,file.path(out,'multiple_imputation_sensitivity_each_imputation.csv'),row.names=FALSE)
write.csv(pooled,file.path(out,'multiple_imputation_sensitivity_pooled.csv'),row.names=FALSE)
write.csv(bind_rows(preds),file.path(out,'multiple_imputation_predictions.csv'),row.names=FALSE)
write.csv(bind_rows(split_rows),file.path(out,'multiple_imputation_sample_partitions.csv'),row.names=FALSE)
write.csv(bind_rows(logs),file.path(out,'multiple_imputation_logged_events.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo_MI.txt'))
