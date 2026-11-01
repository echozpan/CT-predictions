# Install the R packages required by the analysis scripts.
packages <- c('dplyr','ggplot2','glmnet','haven','lme4','mice','naniar',
 'pROC','RANN','ranger','tidyr','MASS')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) install.packages(missing,repos='https://cloud.r-project.org')
