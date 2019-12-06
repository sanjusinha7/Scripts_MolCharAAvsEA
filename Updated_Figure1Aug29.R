# Aug28: Updated Figure 1
##################
#NCIMDProcessing
##################
mat=read.csv('/users/sinhas8/Project_Chromotrypsis/2.Data/Corrected_NCIMD_HRD_by_LOH_and_GI.csv')
mat=mat[mat$hist=='adeno' | mat$hist =='sq',]
mat$hist=as.character(mat$hist)
mat$hist=factor(mat$hist)
levels(mat$hist)=c('LUAD', 'LUSC')
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
scaling_cancerType<-function(quan1=gi, quan2=hist){
  unlist(lapply(split(quan1, quan2), function(x) range01(x)))
}
mat=mat[order(mat$hist),]
mat$Normalized_gi=scaling_cancerType(mat$CNV_burden, mat$hist)
mat$Normalized_hrd.loh=scaling_cancerType(mat$HRD_by_LOH, mat$hist)
mat$Normalized_Chromothripsis_Presence=scaling_cancerType(mat$chtp_Quan, mat$hist)
##################
# Getting Legend
##################
ylim1=c(0,1)
leg=get_legend(ggplot(tcga, aes(x=race,y=Normalized_gi, fill=race))+
                 geom_violin(trim=FALSE)+
                 geom_boxplot(width=0.3, fill="white")+
                 labs(title="GI in PanCan (TCGA)",x="Race", y = "Scaled GI")+
                 theme_classic(base_size = 25)+
                 stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                                    label.y=ylim1[2], size = 10)+
                 theme(legend.box = "horizontal", legend.position = "bottom",
                       legend.text=element_text(size=20), legend.title=element_text(size=25),
                       legend.key.size = unit(2,"line"))+
                 guides(fill=guide_legend(title="Race"))+
                 scale_fill_manual(values = c('#e41a1c','#377eb8'),
                                   labels=c('African Americans', 'European Americans'),
                                   name="Race")
)

##################
# Merging all together
##################
df1=data.frame(Value=mat[mat$hist=='LUSC',]$Normalized_gi, Type='Genomic Instability', Cohort='LUSC NCIMD',
               race=mat[mat$hist=='LUSC',]$race)
df2=data.frame(Value=mat[mat$hist=='LUAD',]$Normalized_gi, Type='Genomic Instability', Cohort='LUAD NCIMD',
               race=mat[mat$hist=='LUAD',]$race)

df3=data.frame(Value=tcga[lung_can,]$Normalized_gi, Type='Genomic Instability', Cohort='LUSC TCGA',
               race=tcga[lung_can,]$race)


df4=data.frame(Value=mat[mat$hist=='LUSC',]$Normalized_hrd.loh, Type='HR-deficiency', Cohort='LUSC NCIMD',
               race=mat[mat$hist=='LUSC',]$race)
df5=data.frame(Value=tcga[lung_can,]$Normalized_hrd.loh, Type='HR-deficiency', Cohort='LUSC TCGA', race=tcga[lung_can,]$race)
df6=data.frame(Value=mat[mat$hist=='LUAD',]$Normalized_hrd.loh, Type='HR-deficiency', Cohort='LUAD NCIMD',
               race=mat[mat$hist=='LUAD',]$race)
mat$chtp_presence=as.numeric(as.logical(mat$chtp_Quan))

df7=data.frame(aggregate(chtp_presence ~ race, function(x) sum(x)/length(x), data = mat[mat$hist=='LUSC',]), Cohort='LUSC NCIMD')
colnames(df7)[2]='CHTP_Canonical_Definition_Presence'

df8=data.frame(aggregate(chtp_presence ~ race, function(x) sum(x)/length(x), data = mat[mat$hist=='LUAD',]), Cohort='LUAD NCIMD')
colnames(df8)[2]='CHTP_Canonical_Definition_Presence'

df9=data.frame(aggregate(CHTP_Canonical_Definition_Presence ~ race, function(x) sum(x)/length(x), data = tcga[lung_can,]), Cohort='LUSC TCGA')
colnames(df9)[2]='CHTP_Canonical_Definition_Presence'
##
df_GI=rbind(df1, df2, df3)
df_HRD=rbind(df4, df5, df6)
df_CHTP=cbind(rbind(df7, df8, df9),Type='Chromothripsis')
##############
# GI
##############
GI_plot=ggplotGrob(
  ggplot(df_GI, aes(x=race,y=Value, fill=race))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.3, fill="white")+
    labs(x="", y = "Scaled GI")+
    theme_bw(base_size = 25)+
    stat_compare_means(method = "wilcox.test", 
                       #label = "p.signif", 
                       label = 'p.format', 
                       label.x = 1.2,
                       label.y = 1.33,
                       #label.y=ylim1[2],
                       size = 8)+
    guides(fill=FALSE)+
    scale_fill_brewer(palette=colorType_Set)+
    facet_grid(Cohort~Type, switch = 'y')+
    theme(#strip.background.y = element_blank(),
      #strip.text.y = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
  #               coord_cartesian(ylim = ylim1*1.15))
)
##############
# HRD
##############
HRD_plot=ggplotGrob(
  ggplot(df_HRD, aes(x=race,y=Value, fill=race))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.3, fill="white")+
    labs(x="Race", y = "Scaled HR-Deficiency")+
    theme_bw(base_size = 25)+
    stat_compare_means(method = "wilcox.test", 
                       #label = "p.signif", 
                       label = 'p',
                       label.x = 1.25,
                       #label.y = 1.3,
                       #label.y=ylim1[2],
                       size = 8)+
    guides(fill=FALSE)+
    scale_fill_brewer(palette=colorType_Set)+
    facet_grid(Cohort~Type)+
    theme(strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
)
##############
# CHTP
##############
###Approprite test for barplot

Pvalues=data.frame(p_values=c(P1=round(fisher.test(t(cbind(table(mat$chtp_presence[mat$hist=='LUSC' & mat$race=='EA']),
                                                         table(mat$chtp_presence[mat$hist=='LUSC' & mat$race=='AA']))),
                                                   alternative = 'g')[1]$p.value, 2),
                              P2=round(fisher.test(cbind(table(mat$chtp_presence[mat$hist=='LUAD' & mat$race=='EA']),
                                                         table(mat$chtp_presence[mat$hist=='LUAD' & mat$race=='AA'])),
                                                   alternative = 'g')[1]$p.value, 2),
                              P3=round(fisher.test(cbind(table(tcga$CHTP_Canonical_Definition_Presence[tcga$hist=='LUSC' & tcga$race=='EA']),
                                                         table(tcga$CHTP_Canonical_Definition_Presence[tcga$hist=='LUSC' & tcga$race=='AA'])),
                                                   alternative = 'g')[1]$p.value, 2)),
                   xstar=c(1.5), ystar=c(0.5, 0.5, 0.55),
                   Cohort=c('LUSC NCIMD', 'LUAD NCIMD', 'LUSC TCGA'),
                   Type=rep('Chromothripsis',3),
                   race=rep('AA', 3)
                   )
#Pvalues$p_values=ifelse(Pvalues$p_values<0.05, '*','ns')

###
df_counts = cbind(rbind(aggregate(X ~ race+hist, mat, function(x) length(x))[c(3, 4, 1, 2),],
                        aggregate(X ~ race+hist, tcga[tcga$hist=='LUSC',], function(x) length(x))),
                  df_CHTP[,-1], x_coor=c(1, 2, 1, 2, 1, 2) )
CHTP_plot=ggplotGrob(
  ggplot(df_CHTP, aes(y=CHTP_Canonical_Definition_Presence, x=race, fill=race))+
    geom_bar(stat = 'identity')+
    labs(x="", y = "Chromothripsis Frequency")+
    theme_bw(base_size = 25)+
    guides(fill=FALSE)+
    scale_fill_brewer(palette=colorType_Set)+
    facet_grid(Cohort~Type, scales = 'free')+
    geom_text(data=Pvalues, aes(x=xstar, y=ystar, label=paste0('p=', p_values)), size=8)+
    geom_text(data=df_counts, aes(label=paste('n=',X, sep=''), x=x_coor), vjust=-0.3, size=7)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background.y = element_blank(),
          strip.text.y = element_blank())
)

##############
# Put all together
##############
##############
# Vector form
##############
setEPS()
postscript('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Figure1_vectorForm.eps',
           height=12, width = 12)
plot_grid(leg,
          plot_grid(GI_plot, HRD_plot, CHTP_plot, align='h', nrow=1,
                    labels = 'AUTO', label_size = 35, rel_widths = c(14/45, 13/45, 10/45) ),
          nrow=2, ncol=1, rel_heights = c(3/50, 47/50))
dev.off()

