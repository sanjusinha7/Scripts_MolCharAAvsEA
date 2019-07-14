####################################
##Aim:: AA specific Novel Drivers
####################################
####################################
##Loading datasets and Pre-Processing
####################################
##Below are frequency files.
require(survival)
err_handle<-function(x){ tryCatch(x, error=function(e){NA}) }
setwd('/Users/sinhas8/')
Freq=read.csv('./Project_Chromotrypsis/2.Data/Novel_AA_specificRecurr_Regions.csv')
Freq_LUAD=read.csv('./Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Freq_LUAD.csv')
Freq_LUSC=read.csv('./Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Freq_LUSC.csv')
###Starting with ~3400 sq and ~2500 candidate recurrent genes
Candidate_AASpecific_Driver_Surv<-function(GeneName, hist, event_type){
  cancer_type_samples=which(mat$hist==hist)
  if(event_type=='Amp'){
    Gene_CNV=CNV[GeneName,na.omit(match(mat$acc[cancer_type_samples], colnames(CNV)))]
    Gene_Event_Dist=ifelse(Gene_CNV>0.1, 'Amp', 'Rest')
  } else{
    Gene_CNV=CNV[GeneName,na.omit(match(mat$acc[cancer_type_samples], colnames(CNV)))]
    Gene_Event_Dist=ifelse(Gene_CNV< -0.1, 'Del', 'Rest')
  }
  
  mat_FORSurvival=mat[!is.na(match(mat$acc, colnames(Gene_CNV))),]
  summary(coxph(Surv(mat_FORSurvival$survival, 
                          mat_FORSurvival$lungcancer_death_all_years)~ factor(unlist(Gene_Event_Dist))))$coefficients[c(1,5)]
  }

Freq=rbind(data.frame(Freq_LUSC, hist='LUSC'),
           data.frame(Freq_LUAD, hist='LUAD'))
colnames(Freq); dim(Freq)
dim(Freq)
df=t(apply(Freq, 1, function(x) Candidate_AASpecific_Driver_Surv(GeneName=as.character(unlist(x[1])),
                                                    hist=as.character(unlist(x[13])),
                                                    event_type=as.character(unlist(x[5])) ))) 
rownames(df)=as.character(Freq$X)
colnames(df)=c('Bcox','pvalue')
df=data.frame(df)
Pot_Novel_Driver_GeneList=rownames(df)
Pot_Novel_Driver_corr=data.frame(do.call(rbind, sapply(as.character(unlist(Pot_Novel_Driver_GeneList)), function(x)
  err_handle(cor.test(unlist(Exp_matched[x,]), unlist(CNV_matched[x,]), method = c('pearson') )[c(3,4)])) ) )
Pot_NDC_df=Freq
Pot_Novel_Driver_corr=Pot_Novel_Driver_corr[which(df$pvalue<0.1 & (Pot_Novel_Driver_corr$p.value<0.1 & Pot_Novel_Driver_corr$estimate>0.1)),]
GOI=as.character(Pot_NDC_df$X)

write.csv(cbind(Freq, df), './Project_Chromotrypsis/2.Data/AA_Spec_Pot_NovelDrivers.csv')
####################################
##Function calling and Plotting
####################################
Candidate_AASpecific_Driver<-function(GeneName, hist, event_type){
  cancer_type_samples=which(mat$hist==hist)
  if(event_type=='Amp'){
    Gene_CNV=CNV[GeneName,na.omit(match(mat$acc[cancer_type_samples], colnames(CNV)))]
    Gene_Event_Dist=ifelse(Gene_CNV>0.1, 'Amp', 'Rest')
  } else{
    Gene_CNV=CNV[GeneName,na.omit(match(mat$acc[cancer_type_samples], colnames(CNV)))]
    Gene_Event_Dist=ifelse(Gene_CNV< -0.1, 'Del', 'Rest')
  }
  
  mat_FORSurvival=mat[!is.na(match(mat$acc, colnames(Gene_CNV))),]
  df=data.frame(mat_FORSurvival, Gene_Event_Dist=factor(unlist(Gene_Event_Dist)))
  fit <- survfit(Surv(survival, 
                      lungcancer_death_all_years) ~ Gene_Event_Dist, 
                 data = df)
  ggsurv<-ggsurvplot(fit, 
                data = df, 
                pval = TRUE,
                xlab = "Time in days",
                ggtheme = theme_bw(),
                title=GeneName,
                legend.labs=levels(df$Gene_Event_Dist),
                )
  ggsurv$plot <- ggsurv$plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  ggsurv
}
temp=apply(Freq[match(GOI, Freq$X),], 
           1, function(x) Candidate_AASpecific_Driver(GeneName=as.character(unlist(x[1])),
                                                                                 hist=as.character(unlist(x[13])),
                                                                                 event_type=as.character(unlist(x[5])) ))
# 
# tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Surv_plot_ofGOI.tif', 
#      width=600, height=600)
# plotA=grid.arrange(grobs=arrange_ggsurvplots(temp, 
#                     nrow = 1,
#                     ncol = length(temp)), 
#                    bottom=textGrob("Time in days", gp=gpar(fontsize=20)),
#                    left=textGrob("Survival Probability", gp=gpar(fontsize=20), rot = 90, vjust = 1))
# dev.off()
####################################
##Aim:: Plot their Expression
####################################
Pot_Novel_Driver_GeneList=GOI
Pot_Novel_Driver_corr=t(sapply(as.character(unlist(Pot_Novel_Driver_GeneList)), function(x)
  cor.test(unlist(Exp_matched[x,]), unlist(CNV_matched[x,]), method = c('pearson') )[c(3,4)]))
Pot_Novel_Driver_corr=Pot_Novel_Driver_corr[!duplicated(rownames(Pot_Novel_Driver_corr)),]
plot_CNV_Exp<-function(GeneName, P_value, Rho, K, legend_position="none"){
  df1=data.frame(Expression=unlist(Exp_matched_scaled[GeneName,]), 
                 CNV=factor(xtile(unlist(CNV_matched[GeneName,]), cutpoints = c(-0.5, -0.1, 0.1, 0.5))) )
  ylim1 = boxplot.stats(df1$Expression)$stats[c(1, 5)]
  xlim1 = boxplot.stats(as.numeric(df1$CNV))$stats[c(1, 5)]
  Figure5=ggplot(data=df1, aes(y=Expression, x=CNV, fill=CNV))+
    geom_boxplot(outlier.colour = NA, outlier.fill = NA)+
    theme_bw(base_size = 15)+
    labs(x='', y='')+
    scale_fill_manual(values = c("Dark green", "green", "grey", 'Pink', 'Dark Red'))+
    theme(legend.position = legend_position, axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    coord_cartesian(ylim = c(ylim1[1]*1.05, ylim1[2]*K))
#    ggtitle(paste(GeneName, ifelse(P_value<0.05, '*', ''), sep='') )
  Figure5
}
##########
###Figure 5B:: Creating Legend
##########
GeneName='TLR6'
legend_position="top"; K=2; P_value=as.numeric(unlist(Pot_Novel_Driver_corr[4,2]))
df1=data.frame(Expression=unlist(Exp_matched_scaled[GeneName,]), 
               CNV=factor(xtile(unlist(CNV_matched[GeneName,]), cutpoints = c(-0.5, -0.1, 0.1, 0.5))) )
ylim1 = boxplot.stats(df1$Expression)$stats[c(1, 5)]
xlim1 = boxplot.stats(as.numeric(df1$CNV))$stats[c(1, 5)]
legend=get_legend(
  ggplot(data=df1, aes(y=Expression, x=CNV, fill=CNV))+
    geom_boxplot(outlier.colour = NA, outlier.fill = NA)+
    #      stat_smooth(method = 'lm', se = FALSE)+
    theme_bw(base_size = 30)+
    labs(x='', y='')+
    scale_fill_manual(values = c("Dark green", "green", "grey", 'Pink', 'Dark Red'),
                      labels=c("Deep Del", "Del", "No Change", 'Amp', 'High Amp'))+
    theme(legend.position = legend_position,
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text = element_text(size=15),
          legend.title = element_blank(),
          legend.key.size = unit(2,"line")
    )+
    guides(fill=guide_legend(nrow=3,byrow=TRUE))+
    coord_cartesian(ylim = c(ylim1[1]*1.05, ylim1[2]*K))+
    ggtitle(paste(GeneName, ifelse(P_value<0.05, '*', ''), sep='') )
)
###############
##Putting them together:: Figure 5B with Legend
###############
Plotexp_forSurvPlot=lapply(rownames(Pot_Novel_Driver_corr), function(x) plot_CNV_Exp(GeneName=x,
                                                                              P_value=Pot_Novel_Driver_corr[x,1],
                                                                              Rho=Pot_Novel_Driver_corr[x,2],
                                                                              K=1))

Comp_Figure_5B=plot_grid(legend,
                         grid.arrange(grobs=Plotexp_forSurvPlot, 
                                      ncol=length(Plotexp_forSurvPlot), nrow=1, 
                                      bottom=textGrob("CNV", gp=gpar(fontsize=20)),
                                      left=textGrob("Expression", gp=gpar(fontsize=20), rot = 90, vjust = 1)
                                      ),
                         ncol=2,
                         rel_widths = c(1/10,9/10)
)
# dev.off()
###############
##Putting them together:: Complete FIgure 5
###############
lay=rbind(c(NA,1,1,1,1,1,1,1,1,1),
          c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2))
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/FigureX_SurvandExp_Correct.tif', 
     width=3000, height=600)
grid.arrange(plotA, 
             Comp_Figure_5B, 
             nrow=2,
             layout_matrix = lay
)
dev.off()

##
Pot_Novel_Drivers=Freq[match(GOI, Freq$X),]
