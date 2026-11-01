# Check that held-out predictors do not affect training imputations.
suppressPackageStartupMessages(library(mice))
args <- commandArgs(TRUE)
code_dir <- if(length(args))args[1] else '.'
expr <- parse(file.path(code_dir,'multiple_imputation.R'))
env <- new.env(parent=globalenv());env$pvs <- paste0('PV',1:5,'CT');env$logs <- list()
for(e in expr) if(is.call(e) && identical(e[[1]],as.name('<-')) &&
  identical(e[[2]],as.name('pmm')))eval(e,env)
set.seed(128)
dat <- data.frame(a=rnorm(250),b=rnorm(250),c=rnorm(250))
dat$b <- dat$a+dat$b
test <- seq_len(250)>200
dat$a[sample(1:200,20)] <- NA
dat$b[sample(201:250,8)] <- NA
altered <- dat;altered$c[test] <- altered$c[test]+100
a <- env$pmm(dat,test,2,5678)
b <- env$pmm(altered,test,2,5678)
stopifnot(all(vapply(1:2,function(i)isTRUE(all.equal(a[[i]][!test,],b[[i]][!test,])),logical(1))))
cat('PASS: altering held-out observed predictors leaves training imputations identical.\n')
