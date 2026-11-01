# Extract the public-use variables required for CT prediction analyses.
suppressPackageStartupMessages({library(dplyr); library(haven)})
args <- commandArgs(TRUE)
if (length(args) != 2) stop('Supply the official SPSS data directory and the extract destination.')
raw_dir <- args[1]
destination <- args[2]
dir.create(file.path(destination, 'raw'), recursive=TRUE, showWarnings=FALSE)
codes <- c('DEU','DNK','FIN','FRA','KOR','LUX','PRT','USA')
student <- c('IDCNTRY','IDSCHOOL','IDSTUD','TOTWGTS','S_SEX','S_AGE','S_HISEI',
 'S_HOMLIT','S_IMMIG','S_SPECEFF','S_GENEFF','S_USECOM','S_USEINF','S_ACCONT',
 'S_SPECLASS','S_GENCLASS','S_CODLRN','S_ICTFUT','S_ICTLRN','S_ICTNEG',
 'S_ICTPOS','S_GENACT','S_USESTD','S_NISB',paste0('PV',1:5,'CIL'),paste0('PV',1:5,'CT'))
school <- c('IDSCHOOL','P_RATTCH','P_EXPLRN','P_PRIORH','P_PRIORS','P_VWICT',
 'P_ICTUSE','P_ICTCOM','P_EXPTCH')
merged <- dictionary <- list()
for (code in codes) {
 sf <- paste0('BSG',code,'I2.sav'); pf <- paste0('BCG',code,'I2.sav')
 s <- read_sav(file.path(raw_dir,sf),col_select=all_of(student))
 p <- read_sav(file.path(raw_dir,pf),col_select=all_of(school))
 stopifnot(!anyDuplicated(p$IDSCHOOL))
 attr(s,'label') <- 'ICILS 2018 public-use student analysis extract'
 attr(p,'label') <- 'ICILS 2018 public-use school analysis extract'
 write_sav(s,file.path(destination,'raw',sf))
 write_sav(p,file.path(destination,'raw',pf))
 d <- left_join(s,p,by='IDSCHOOL')
 stopifnot(nrow(d)==nrow(s))
 dictionary[[code]] <- data.frame(country_code=code,variable=names(d),
  label=vapply(d,function(x){z<-attr(x,'label');if(is.null(z))'' else z},character(1)),
  level=ifelse(names(d) %in% setdiff(school,'IDSCHOOL'),'School','Student'),row.names=NULL)
 d[] <- lapply(d, function(x) as.numeric(zap_labels(x)))
 d$COUNTRY_CODE <- code
 d$SCHOOL_UID <- paste(code,d$IDSCHOOL,sep='_')
 d$ROW_UID <- paste(code,d$IDSTUD,sep='_')
 predictors <- c(setdiff(student,c('IDCNTRY','IDSCHOOL','IDSTUD','TOTWGTS',
  'S_SEX','S_HISEI','S_HOMLIT','S_IMMIG',paste0('PV',2:5,'CIL'),paste0('PV',1:5,'CT'))),setdiff(school,'IDSCHOOL'))
 d$complete_case <- complete.cases(d[,c(predictors,paste0('PV',1:5,'CT'))])
 merged[[code]] <- d
}
con <- gzfile(file.path(destination,'analysis_data.csv.gz'),'wt')
write.csv(bind_rows(merged),con,row.names=FALSE,na='')
close(con)
write.csv(bind_rows(dictionary),file.path(destination,'source_variable_labels.csv'),row.names=FALSE,na='')
cat('Extracted',sum(vapply(merged,nrow,integer(1))),'student records across',length(codes),'national samples.\n')
