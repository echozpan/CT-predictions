# Download the official ICILS 2018 public-use SPSS database.
args <- commandArgs(TRUE)
destination <- if (length(args)) args[1] else 'data/official'
dir.create(destination,recursive=TRUE,showWarnings=FALSE)
options(timeout=max(600,getOption('timeout')))
url <- paste0('https://www.iea.nl/sites/default/files/data-repository/ICILS/',
 'ICILS2018/ICILS2018_IDB_SPSS.zip')
archive <- file.path(destination,'ICILS2018_IDB_SPSS.zip')
message('Source and terms: https://www.iea.nl/data-tools/repository/icils')
download.file(url,archive,mode='wb')
unzip(archive,exdir=destination)
