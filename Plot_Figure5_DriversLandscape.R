########################
##Figure 5A:: Distribution of COSMIC genes Alteration Freq
########################
library(grid)
library(gridExtra)
require(ggrepel)
require(cowplot)
df_Freq=cbind(Amp=read.csv('/Users/sinhas8/Project_Chromotrypsis/4.Results/GISTIC_results/Complete_NCIMD/Gen_Level_Freq_COSMIC/gene_Amp_Freq_COSMIC.csv'),
              Del=read.csv('/Users/sinhas8/Project_Chromotrypsis/4.Results/GISTIC_results/Complete_NCIMD/Gen_Level_Freq_COSMIC/gene_Del_Freq_COSMIC.csv')[,-1])
rownames(df_Freq) = df_Freq[,1];df_Freq=df_Freq[,-1] 
cosmic=read.csv('/Users/sinhas8/Project_Chromotrypsis/2.Data/COSMIC.csv', sep='\t')
lung_cosmic=cosmic[sapply(seq(nrow(cosmic)), function(x) 
  sum(((grepl('lung', sapply(c(cosmic[x,]), as.character), ignore.case = TRUE)|
          grepl('NSCLC', sapply(c(cosmic[x,]), as.character), , ignore.case = TRUE))|
         grepl('multiple', sapply(c(cosmic[x,]), as.character), , ignore.case = TRUE) )|
        grepl('JAK2', sapply(c(cosmic[x,]), as.character), , ignore.case = TRUE) )>0),]
lung_cosmic=rbind(lung_cosmic, cosmic[cosmic$Gene.Symbol=='PTEN',])
lung_cosmic$Gene_Type=NA
lung_cosmic$Gene_Type[grep('TSG',lung_cosmic$Role.in.Cancer, ignore.case = T)]='TSG'
lung_cosmic$Gene_Type[grep('oncogene',as.character(unlist(lung_cosmic$Role.in.Cancer)), ignore.case = T)]='Oncogene'
lung_cosmic$Gene_Type[grep('oncogene, TSG',
                           as.character(unlist(lung_cosmic$Role.in.Cancer)), ignore.case = T)]='Oncogene and TSG'
# dim(lung_cosmic); dim(df_Freq)
df_Freq_lung=df_Freq[match(lung_cosmic$Gene.Symbol, rownames(df_Freq) ),]
df_Freq_lung$Gene_Type=lung_cosmic$Gene_Type
df_Freq_lung=cbind(rownames(df_Freq_lung), df_Freq_lung)

Amp_df_Freq_lung=df_Freq_lung[grep('Oncogene',df_Freq_lung$Gene_Type),1:5]
Del_df_Freq_lung=df_Freq_lung[grep('TSG',df_Freq_lung$Gene_Type),c(1,6:9)]

Amp_df_Freq_lung_LUAD=data.frame(Amp_df_Freq_lung[,c(1,2,3)], Histology='LUAD', Event_Type='Amp')
Amp_df_Freq_lung_LUSC=data.frame(Amp_df_Freq_lung[,c(1,4,5)], Histology='LUSC', Event_Type='Amp')

Del_df_Freq_lung_LUAD=data.frame(Del_df_Freq_lung[,c(1,2,3)], Histology='LUAD', Event_Type='Del')
Del_df_Freq_lung_LUSC=data.frame(Del_df_Freq_lung[,c(1,4,5)], Histology='LUSC', Event_Type='Del')
colnames(Amp_df_Freq_lung_LUAD)[1:3]=c('GeneName','AA_Freq','EA_Freq')
colnames(Del_df_Freq_lung_LUSC)[1:3]=c('GeneName','AA_Freq','EA_Freq')
colnames(Amp_df_Freq_lung_LUSC)[1:3]=c('GeneName','AA_Freq','EA_Freq')
colnames(Del_df_Freq_lung_LUAD)[1:3]=c('GeneName','AA_Freq','EA_Freq')
Freq_lung=rbind(Amp_df_Freq_lung_LUAD, Amp_df_Freq_lung_LUSC, 
                Del_df_Freq_lung_LUAD, Del_df_Freq_lung_LUSC)
################
###Finding Recurrence
################
Recurrence_Type<-function(GeneName='CDKN2A', hist='LUSC', Event_Type){
  #For AA
  Race='AA'
  df1=read.csv(paste('/Users/sinhas8/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/',
                     hist, '_',Race,'/', 'freq_of_cytobands_wd_Genes.csv', sep=''))
  ReccurntIn_AA=sum(grepl(GeneName,df1$Genes[grep(Event_Type,df1$Unique.Name)]))>0
  #For EA
  Race='EA'
  df1=read.csv(paste('/Users/sinhas8/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/',
                     hist, '_',Race,'/', 'freq_of_cytobands_wd_Genes.csv', sep=''))
  ReccurntIn_EA=sum(grepl(GeneName,df1$Genes[grep(Event_Type,df1$Unique.Name)]))>0
  if(ReccurntIn_AA & ReccurntIn_EA){ return('Recurrent in Both')
  } else if(ReccurntIn_AA) { return('Recurrent in AA Only')
  } else if(ReccurntIn_EA) { return('Recurrent in EA Only') 
  } else { return('Recurrent in None') } 
}

Freq_lung$Significance=apply(Freq_lung, 1, function(x) 
  err_handle(Recurrence_Type(GeneName=as.character(unlist(x[1])),
                  hist=as.character(unlist(x[4])),
                  Event_Type=as.character(unlist(x[5])) ) ))
################
###Plotting
################
Freq_lung_wdoutNA=na.omit(Freq_lung)
Figure5A=ggplotGrob(
  ggplot(Freq_lung_wdoutNA, aes(x=AA_Freq, y=EA_Freq, label = GeneName, color=Significance))+
    geom_point(position=position_jitter(h=0.01,w=0.01), size=2.5)+ 
    geom_abline(intercept = 0, slope = 1, color='blue')+
    geom_text_repel(size=8, show.legend = FALSE)+ 
    facet_grid(Histology~Event_Type)+
    theme_bw(base_size = 24)+
    labs(x='Frequency in AA', y='Frequency in EA', color="SCNA recurrence significance")+
    #ggtitle('Lung cancer Drivers Alteration Freq across race')+
    scale_color_manual(values=c("blue", "orange", "red", "black"))+ 
    guides(color=guide_legend(nrow=2, byrow=TRUE, override.aes = list(size = 3) ))+
    #   guides(color=guide_legend()+
    theme(legend.position="top", legend.text=element_text(size=22),
          legend.title=element_text(size=24),
          axis.title = element_text(size=28))
)
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Figure5_test1_LUSC_May7.tif',
     height = 1200, width = 1200)
plot(Figure5A)
dev.off()
##########
###Figure 5B:: Showing Exp and CNV Changes
##########
source('../3.Tools/define_Exp_matched_Scale.R')
Exp_matched_scaled=t(Exp_matched_scaled)
CNV_matched_scaled=t(CNV_matched_scaled)
Drivers=Freq_lung_wdoutNA[-grep('None',Freq_lung_wdoutNA$Significance),]
Driver_GeneList=unique(Drivers$GeneName)

Driver_corr=t(sapply(as.character(unlist(Driver_GeneList)), function(x)
  cor.test(unlist(Exp_matched[x,]), unlist(CNV_matched[x,]), method = c('pearson') )[c(3,4)]))
Driver_corr=Driver_corr[!duplicated(rownames(Driver_corr)),]
plot_CNV_Exp<-function(GeneName, P_value, Rho, K, legend_position="none"){
  df1=data.frame(Expression=unlist(Exp_matched_scaled[GeneName,]), 
                 CNV=factor(xtile(unlist(CNV_matched[GeneName,]), cutpoints = c(-0.5, -0.1, 0.1, 0.5))) )
  ylim1 = boxplot.stats(df1$Expression)$stats[c(1, 5)]
  xlim1 = boxplot.stats(as.numeric(df1$CNV))$stats[c(1, 5)]
  Figure5=ggplot(data=df1, aes(y=Expression, x=CNV, fill=CNV))+
    geom_boxplot(outlier.colour = NA, outlier.fill = NA)+
    #      stat_smooth(method = 'lm', se = FALSE)+
    theme_bw(base_size = 20)+
    labs(x='', y='')+
    scale_fill_manual(values = c("Dark green", "green", "grey", 'Pink', 'Dark Red'))+
    theme(legend.position = legend_position, axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    coord_cartesian(ylim = c(ylim1[1]*1.05, ylim1[2]*K))+
    ggtitle(paste(GeneName, ifelse(P_value<0.05, '*', ''), sep='') )
  Figure5
}
##########
###Figure 5B:: Creating Legend
##########
GeneName='CDKN2A'
legend_position="top"; K=2; P_value=0.1
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
          legend.text = element_text(size=20),
          legend.title = element_blank(),
          legend.key.size = unit(2,"line")
    )+
    guides(fill=guide_legend(nrow=3,byrow=TRUE))+
    coord_cartesian(ylim = c(ylim1[1]*1.05, ylim1[2]*K))+
    ggtitle(paste(GeneName, ifelse(P_value<0.05, '*', ''), sep='') )
)
###############
##Putting them together:: Supp Figure of Expression
###############
Supp_Figure5B_plist=c(lapply(rownames(Driver_corr)[c(2, 4, 5, 7, 13)], function(x) plot_CNV_Exp(GeneName=x,
                                                                                           P_value=Driver_corr[x,1],
                                                                                           Rho=Driver_corr[x,2],
                                                                                           K=4.5)),
                 lapply(rownames(Driver_corr)[-c(2, 4, 5, 7, 13)], function(x) plot_CNV_Exp(GeneName=x,
                                                                                            P_value=Driver_corr[x,1],
                                                                                            Rho=Driver_corr[x,2],
                                                                                            K=1.5)))
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Supp_Figure_Expressioncorr_Sep1.tif', 
     width=1500, height=1200)
plot_grid(legend,
          do.call("grid.arrange", c(Supp_Figure5B_plist, ncol=6, nrow=3)),
          nrow=2, rel_heights = c(1/10, 9/10))
dev.off()
###############
##Putting them together:: Figure 5B with Legend
###############
Exp_matched_scaled['CDKN1A',]
Figure5B_plist=c(lapply(rownames(Driver_corr)[c(4)], function(x) plot_CNV_Exp(GeneName=x,
                                                                              P_value=Driver_corr[x,1],
                                                                              Rho=Driver_corr[x,2],
                                                                              K=4.5)),
                 lapply(rownames(Driver_corr)[c(7)], function(x) plot_CNV_Exp(GeneName=x,
                                                                              P_value=Driver_corr[x,1],
                                                                              Rho=Driver_corr[x,2],
                                                                              K=4.5)),
                 lapply(rownames(Driver_corr)[c(18)], function(x) plot_CNV_Exp(GeneName=x,
                                                                              P_value=Driver_corr[x,1],
                                                                              Rho=Driver_corr[x,2],
                                                                              K=1)))

# tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Figure5B.tif', 
#        width=1200, height=1200)
grobs <- list()
widths <- list()
for (i in 1:length(Figure5B_plist)){
  grobs[[i]] <- ggplotGrob(Figure5B_plist[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}
maxwidth <- do.call(grid::unit.pmax, widths)
for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}
Comp_Figure_5B=plot_grid(legend,
          grid.arrange(grobs=grobs, ncol=1, nrow=3, 
                       bottom=textGrob("CNV", gp=gpar(fontsize=25)),
                       left=textGrob("Expression", gp=gpar(fontsize=25), rot = 90, vjust = 1)),
          nrow=2,
          rel_heights = c(2/10,8/10)
)
# dev.off()
###############
##Putting them together:: Complete FIgure 5
###############
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Comp_Figure5_Sep1.tif', 
     width=1800, height=1200)
plot_grid(Figure5A, 
          plot_grid(NULL, Comp_Figure_5B, NULL, rel_heights = c(1/10, 3/5, 1/10), nrow=3), 
          ncol=2,
          labels='AUTO',  label_size = 35,
          #layout_matrix = lay, 
          rel_widths = c(4/5, 1/5)
)
dev.off()

###############
##Adding Survival
###############
require(survival)
Driver_OI=names(which(unlist(Driver_corr[,1])<0.1))
CNV_FORSurvival=CNV[,na.omit(match(mat$acc[mat$hist=='LUSC'], colnames(CNV)))]
mat_FORSurvival=mat[!is.na(match(mat$acc, colnames(CNV_FORSurvival))),]
temp=apply(CNV_FORSurvival[Driver_OI,], 1, function(x) 
  summary(coxph(Surv(mat_FORSurvival$survival, mat_FORSurvival$lungcancer_death_5_years)~ x))$coefficients[c(1,5)])
temp
CNV_FORSurvival=CNV[,na.omit(match(mat$acc[mat$hist=='LUAD'], colnames(CNV)))]
mat_FORSurvival=mat[!is.na(match(mat$acc, colnames(CNV_FORSurvival))),]
temp=apply(CNV_FORSurvival[Driver_OI,], 1, function(x) 
  summary(coxph(Surv(mat_FORSurvival$survival, mat_FORSurvival$lungcancer_death_5_years)~ x))$coefficients[c(1,5)])
temp

