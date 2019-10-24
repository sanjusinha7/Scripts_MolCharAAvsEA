##Libraries required
require(ggplot2)
require(ggrepel)
setwd('/Users/sinhas8/Project_Chromotrypsis/4.Results/')
##Plot of arm-level events frequemcy distribution
########################
###Figure 2A #in NCIMD
########################
adeno_AA=read.csv('./GISTIC_results/Trimmed_GISTIC_Results/LUAD_AA/broad_significance_results.txt', sep='\t')
adeno_EA=read.csv('./GISTIC_results/Trimmed_GISTIC_Results/LUAD_EA/broad_significance_results.txt', sep='\t')
sq_AA=read.csv('./GISTIC_results/Trimmed_GISTIC_Results/LUSC_EA/broad_significance_results.txt', sep='\t')
sq_EA=read.csv('./GISTIC_results/Trimmed_GISTIC_Results/LUSC_EA/broad_significance_results.txt', sep='\t')

##
df_sq_Amp=data.frame(Arm=sq_AA$Arm, AA_Freq=sq_AA$Amp.frequency, EA_Freq=sq_EA$Amp.frequency, AA_qvalue=sq_AA$Amp.q.value, EA_qvalue=sq_EA$Amp.q.value,
                        CNV_Type='Amp')
df_sq_Del=data.frame(Arm=sq_AA$Arm, AA_Freq=sq_AA$Del.frequency, EA_Freq=sq_EA$Del.frequency, AA_qvalue=sq_AA$Del.q.value, EA_qvalue=sq_EA$Del.q.value,
                        CNV_Type='Del')
##
df_adeno_Amp=data.frame(Arm=adeno_AA$Arm, AA_Freq=adeno_AA$Amp.frequency, EA_Freq=adeno_EA$Amp.frequency, AA_qvalue=adeno_AA$Amp.q.value, EA_qvalue=adeno_EA$Amp.q.value,
                        CNV_Type='Amp')
df_adeno_Del=data.frame(Arm=adeno_AA$Arm, AA_Freq=adeno_AA$Del.frequency, EA_Freq=adeno_EA$Del.frequency, AA_qvalue=adeno_AA$Del.q.value, EA_qvalue=adeno_EA$Del.q.value,
                        CNV_Type='Del')
###Adding recurrence infomation for Adeno
df_adeno_Amp$Sig = 'Recurrent in None'
df_adeno_Amp$Sig[df_adeno_Amp$AA_qvalue<0.05 & df_adeno_Amp$EA_qvalue<0.05] = 'Recurrent in Both'
df_adeno_Amp$Sig[df_adeno_Amp$AA_qvalue<0.05 & df_adeno_Amp$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_adeno_Amp$Sig[df_adeno_Amp$AA_qvalue>0.05 & df_adeno_Amp$EA_qvalue<0.05] = 'Recurrent in EA Only'
df_adeno_Del$Sig = 'Recurrent in None'
df_adeno_Del$Sig[df_adeno_Del$AA_qvalue<0.05 & df_adeno_Del$EA_qvalue<0.05] = 'Recurrent in Both'
df_adeno_Del$Sig[df_adeno_Del$AA_qvalue<0.05 & df_adeno_Del$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_adeno_Del$Sig[df_adeno_Del$AA_qvalue>0.05 & df_adeno_Del$EA_qvalue<0.05] = 'Recurrent in EA Only'

df_adeno=rbind(df_adeno_Amp, df_adeno_Del)
###Adding recurrence infomation for sq
df_sq_Amp$Sig = 'Recurrent in None'
df_sq_Amp$Sig[df_sq_Amp$AA_qvalue<0.05 & df_sq_Amp$EA_qvalue<0.05] = 'Recurrent in Both'
df_sq_Amp$Sig[df_sq_Amp$AA_qvalue<0.05 & df_sq_Amp$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_sq_Amp$Sig[df_sq_Amp$AA_qvalue>0.05 & df_sq_Amp$EA_qvalue<0.05] = 'Recurrent in EA Only'
df_sq_Del$Sig = 'Recurrent in None'
df_sq_Del$Sig[df_sq_Del$AA_qvalue<0.05 & df_sq_Del$EA_qvalue<0.05] = 'Recurrent in Both'
df_sq_Del$Sig[df_sq_Del$AA_qvalue<0.05 & df_sq_Del$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_sq_Del$Sig[df_sq_Del$AA_qvalue>0.05 & df_sq_Del$EA_qvalue<0.05] = 'Recurrent in EA Only'

df_sq=rbind(df_sq_Amp, df_sq_Del)

#Formatting the dfs
df=rbind(cbind(df_adeno, Hist='LUAD'), cbind(df_sq, Hist='LUSC'))
colnames(df)[7]='Significance'
df_NCIMD=df
###Plot of Arm level events
########################
###Figure 2B #in TCGA
########################
adeno_AA=read.csv('gistic_TCGA/LUAD_AA/broad_significance_results.txt', sep='\t')
adeno_EA=read.csv('gistic_TCGA/LUAD_EA/broad_significance_results.txt', sep='\t')
sq_AA=read.csv('gistic_TCGA/LUSC_AA/broad_significance_results.txt', sep='\t')
sq_EA=read.csv('gistic_TCGA/LUSC_EA/broad_significance_results.txt', sep='\t')

##
df_sq_Amp=data.frame(Arm=sq_AA$Arm, AA_Freq=sq_AA$Amp.frequency, EA_Freq=sq_EA$Amp.frequency, AA_qvalue=sq_AA$Amp.q.value, EA_qvalue=sq_EA$Amp.q.value,
                     CNV_Type='Amp')
df_sq_Del=data.frame(Arm=sq_AA$Arm, AA_Freq=sq_AA$Del.frequency, EA_Freq=sq_EA$Del.frequency, AA_qvalue=sq_AA$Del.q.value, EA_qvalue=sq_EA$Del.q.value,
                     CNV_Type='Del')
##
df_adeno_Amp=data.frame(Arm=adeno_AA$Arm, AA_Freq=adeno_AA$Amp.frequency, EA_Freq=adeno_EA$Amp.frequency, AA_qvalue=adeno_AA$Amp.q.value, EA_qvalue=adeno_EA$Amp.q.value,
                        CNV_Type='Amp')
df_adeno_Del=data.frame(Arm=adeno_AA$Arm, AA_Freq=adeno_AA$Del.frequency, EA_Freq=adeno_EA$Del.frequency, AA_qvalue=adeno_AA$Del.q.value, EA_qvalue=adeno_EA$Del.q.value,
                        CNV_Type='Del')
###Adding recurrence infomation for Adeno
df_adeno_Amp$Sig = 'Recurrent in None'
df_adeno_Amp$Sig[df_adeno_Amp$AA_qvalue<0.05 & df_adeno_Amp$EA_qvalue<0.05] = 'Recurrent in Both'
df_adeno_Amp$Sig[df_adeno_Amp$AA_qvalue<0.05 & df_adeno_Amp$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_adeno_Amp$Sig[df_adeno_Amp$AA_qvalue>0.05 & df_adeno_Amp$EA_qvalue<0.05] = 'Recurrent in EA Only'
df_adeno_Del$Sig = 'Recurrent in None'
df_adeno_Del$Sig[df_adeno_Del$AA_qvalue<0.05 & df_adeno_Del$EA_qvalue<0.05] = 'Recurrent in Both'
df_adeno_Del$Sig[df_adeno_Del$AA_qvalue<0.05 & df_adeno_Del$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_adeno_Del$Sig[df_adeno_Del$AA_qvalue>0.05 & df_adeno_Del$EA_qvalue<0.05] = 'Recurrent in EA Only'

df_adeno=rbind(df_adeno_Amp, df_adeno_Del)
###Adding recurrence infomation for sq
df_sq_Amp$Sig = 'Recurrent in None'
df_sq_Amp$Sig[df_sq_Amp$AA_qvalue<0.05 & df_sq_Amp$EA_qvalue<0.05] = 'Recurrent in Both'
df_sq_Amp$Sig[df_sq_Amp$AA_qvalue<0.05 & df_sq_Amp$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_sq_Amp$Sig[df_sq_Amp$AA_qvalue>0.05 & df_sq_Amp$EA_qvalue<0.05] = 'Recurrent in EA Only'
df_sq_Del$Sig = 'Recurrent in None'
df_sq_Del$Sig[df_sq_Del$AA_qvalue<0.05 & df_sq_Del$EA_qvalue<0.05] = 'Recurrent in Both'
df_sq_Del$Sig[df_sq_Del$AA_qvalue<0.05 & df_sq_Del$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_sq_Del$Sig[df_sq_Del$AA_qvalue>0.05 & df_sq_Del$EA_qvalue<0.05] = 'Recurrent in EA Only'

df_sq=rbind(df_sq_Amp, df_sq_Del)

#Formatting the dfs
df=rbind(cbind(df_adeno, Hist='LUAD'), cbind(df_sq, Hist='LUSC'))
colnames(df)[7]='Significance'
df_TCGA=df
df=rbind(data.frame(df_NCIMD, Cohort='NCIMD'),
         data.frame(df_TCGA, Cohort='TCGA') )

########################
#Plot them together
########################
AX=ggplotGrob(
  ggplot(df, aes(x=AA_Freq, y=EA_Freq, color=Significance, label = Arm))+
    geom_point(position=position_jitter(h=0.01,w=0.01), size=2.5)+ 
    geom_abline(intercept = 0, slope = 1, color='blue', linetype='dashed')+
    geom_text_repel(size=8, show.legend = FALSE)+ 
    facet_grid(Cohort ~ Hist+CNV_Type, scales = 'fixed')+
    theme_bw(base_size = 24)+
    scale_color_manual(values=c("red", "orange", "blue", "black"))+ 
    guides(color=guide_legend(nrow=2, byrow=TRUE, override.aes = list(size = 3) ))+
    #   guides(color=guide_legend()+
    theme(legend.position="top", legend.text=element_text(size=20),
          legend.title=element_text(size=22))+
    #lims(x=c(0, max(max(df_adeno$AA_Freq), max(df_adeno$EA_Freq))), y=c(0,max(max(df_adeno$AA_Freq), max(df_adeno$EA_Freq))))	+
    labs(x='Frequency in AA', y='Frequency in EA', color="SCNA recurrence significance")+
    coord_cartesian(xlim = c(0, 0.8), ylim = c(0, 0.8)) 
)

tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/31stAugFigure3_TCGAandNCIMDv1.tif',
     width = 1800, height = 1000)
plot(AX)
dev.off()
########################
#Figure 2B
########################
require(VennDiagram); install.packages('metafolio')
B1=draw.pairwise.venn(area1= 30,
                      area2= 27,
                      cross.area= 12,
    #                  category= c("AA", "EA"),
                      fill= c('#e41a1c','#377eb8'),
                      main='Recurrent',
                      cex= 2,
                      scaled=FALSE,
                      alpha=0.8,
                      cat.cex = 2
)
B2= draw.pairwise.venn(area1= 45,
                       area2= 54,
                       cross.area= 24,
   #                    category= c("AA", "EA"),
                       fill= c('#e41a1c','#377eb8'),
                       main='Recurrent',
                       cex= 2,
                       scaled=FALSE,
                       alpha=0.8,
                       cat.cex = 2, 
                       rotation.degree=180
)
B3= draw.pairwise.venn(area1= 18,
                       area2= 34,
                       cross.area= 10,
                       #category= c("AA", "EA"),
                       fill= c('#e41a1c','#377eb8'),
                       main='Recurrent',
                       cex= 2,
                       scaled=FALSE,
                       alpha=0.8,
                       cat.cex = 2,
                       rotation.degree=180
)
B4= draw.pairwise.venn(area1= 32,
                       area2= 34,
                       cross.area= 12,
                       #category= c("AA", "EA"),
                       fill= c('#e41a1c','#377eb8'),
                       main='Recurrent',
                       cex= 2,
                       scaled=FALSE,
                       alpha=0.8,
                       cat.cex = 2,
                       rotation.degree=180
)
require(cowplot); require(gridExtra)
Figure3_LUAD = plot_grid(grid.arrange(gTree(children=B1), top=textGrob("Amp", gp=gpar(fontsize=20,font=8))),
              grid.arrange(gTree(children=B2), top=textGrob("Del", gp=gpar(fontsize=20,font=8))),
#              grid.arrange(gTree(children=B3), top=textGrob("Recurrent Focal Amp in LUSC", gp=gpar(fontsize=20,font=8))),
 #             grid.arrange(gTree(children=B4), top=textGrob("Recurrent Focal Del in LUSC", gp=gpar(fontsize=20,font=8))),
              nrow=1)
Figure3_LUSC = plot_grid(grid.arrange(gTree(children=B3), top=textGrob("Amp", gp=gpar(fontsize=20,font=8))),
                         grid.arrange(gTree(children=B4), top=textGrob("Del", gp=gpar(fontsize=20,font=8))),
                         nrow=1)
getwd()
tiff('prep_final_figures/Fig3LUAD.tif', width=300, height = 150)
plot(Figure3_LUAD)
dev.off()
tiff('prep_final_figures/Fig3LUSC.tif', width=300, height = 150)
plot(Figure3_LUSC)
dev.off()


########################
#Figure 2C
########################
df1_pp=read.csv('/Users/sinhas8/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Corr_with_Survival_and_Expression_sq.csv')
df1_pp = df1_pp[,-c(2:3, 6:7)]
###
t1=df1_pp[grepl('Del',df1_pp[,2]) | grepl('Del',df1_pp[,3]),c(1,6,7)]
colnames(t1)=c('GeneName','Freq_EA','Freq_AA')
t2=df1_pp[!(grepl('Del',df1_pp[,2]) | grepl('Del',df1_pp[,3])),c(1,4,5)]
colnames(t2)=c('GeneName','Freq_EA','Freq_AA')
df2plot=rbind(t1, t2)
rownames(df2plot)=df2plot[,1]
df2plot=df2plot[,-1]
df2plot_num=apply(df2plot, 2, as.numeric)
rownames(df2plot_num)=rownames(df2plot)
df2plot$Type=c(rep('Del', nrow(df2plot)-2), 
                          rep('Amp', 2))
require('NMF')
Panel3=ggplotGrob(
  ggplot(df2plot, aes(x=Freq_AA, y=Freq_EA, label = rownames(df2plot)))+
    geom_point(position=position_jitter(h=0.01,w=0.01), size=2.5)+ 
    geom_abline(intercept = 0, slope = 1, color='blue')+
    geom_text_repel(size=8, show.legend = FALSE)+ 
    facet_grid(~Type)+
    theme_bw(base_size = 24)+
    labs(x='Frequency in AA', y='Frequency in EA')+
    ggtitle('Frequency of recurrent genes associated with Expression and Survival')
)

##########
###Final Completed Figure
##########
