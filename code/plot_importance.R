# Plot held-out permutation importance and predictor missingness.
suppressPackageStartupMessages({library(dplyr);library(ggplot2)})
out <- commandArgs(TRUE)[1]
imp <- read.csv(file.path(out,'permutation_importance_all_models.csv'))
for(context in unique(imp$context))for(code in c('KOR','DNK')) {
  z <- imp |> filter(.data$context==.env$context,target_code==code) |>
    group_by(model) |> slice_max(delta_MSE,n=6,with_ties=FALSE) |> ungroup()
  z$key <- paste(z$variable,z$model,sep='___')
  z$key <- factor(z$key,levels=z$key[order(z$delta_MSE)])
  g <- ggplot(z,aes(delta_MSE,key,fill=model))+geom_col(show.legend=FALSE)+
    facet_wrap(~model,ncol=2,scales='free_y')+
    scale_y_discrete(labels=function(x)sub('___.*','',x))+
    labs(x='Increase in held-out MSE after permutation',y=NULL,
      title=paste(context,ifelse(code=='KOR','Korea','Denmark'),sep=' - '),
      subtitle='PV1CT; six highest variables per model; three permutations')+
    theme_minimal(base_size=11)+theme(panel.grid.major.y=element_blank())
  stem <- ifelse(context=='Within-country','within','external')
  ggsave(file.path(out,'figures',paste0('importance_',stem,'_',code,'.png')),
    g,width=7,height=7,dpi=300)
}
mr <- read.csv(file.path(out,'missing_rates_by_country_variable.csv'))
mr$variable <- factor(mr$variable,levels=rev(unique(mr$variable)))
g <- ggplot(mr,aes(country,variable,fill=percent_missing))+
  geom_tile(color='white',linewidth=.15)+scale_fill_viridis_c(option='C',direction=-1)+
  labs(x=NULL,y=NULL,fill='% missing',title='Missingness before complete-case selection')+
  theme_minimal(base_size=12)+theme(axis.text.x=element_text(angle=35,hjust=1))
ggsave(file.path(out,'figures','missingness_heatmap.png'),g,width=9,height=8,dpi=300)

p <- read.csv(file.path(out,'model_metrics_pooled_across_plausible_values.csv'))
g <- p |> filter(feature_set=='All predictors') |>
 ggplot(aes(model,RMSE_mean,color=context,group=context))+
 geom_point(size=2.4,position=position_dodge(.25))+
 geom_line(position=position_dodge(.25))+facet_wrap(~target_country)+
 coord_cartesian(ylim=c(0,NA))+labs(x=NULL,y='RMSE (ICILS CT scale points)',color=NULL,
 title='Within-country and external-validation error',
 subtitle='Means across five plausible values; school-blocked target holdouts')+
 theme_minimal(base_size=11)+theme(axis.text.x=element_text(angle=30,hjust=1),legend.position='bottom')
ggsave(file.path(out,'figures','performance_comparison.png'),g,width=8.5,height=5.2,dpi=300)
top <- imp |> group_by(context,target_country,model) |> slice_max(delta_MSE,n=8,with_ties=FALSE) |>
 ungroup() |> mutate(variable=reorder(variable,delta_MSE))
g <- ggplot(top,aes(delta_MSE,variable,fill=model))+geom_col(show.legend=FALSE)+
 facet_grid(target_country~context+model,scales='free_y')+
 labs(x='Increase in MSE after permutation',y=NULL,title='Permutation importance across models and targets')+
 theme_minimal(base_size=8)
ggsave(file.path(out,'figures','feature_importance.png'),g,width=15,height=8,dpi=300)
