# Plot of arm-level events frequemcy distribution
# Libraries required
require(ggplot2)
require(ggrepel)
setwd('/Users/sinhas8/Project_Chromotrypsis/4.Results/')
########################
# Figure 2A #in NCIMD
########################
adeno_AA=read.csv('./GISTIC_results/Trimmed_GISTIC_Results/LUAD_AA/broad_significance_results.txt', sep='\t')
adeno_EA=read.csv('./GISTIC_results/Trimmed_GISTIC_Results/LUAD_EA/broad_significance_results.txt', sep='\t')
sq_AA=read.csv('./GISTIC_results/Trimmed_GISTIC_Results/LUSC_EA/broad_significance_results.txt', sep='\t')
sq_EA=read.csv('./GISTIC_results/Trimmed_GISTIC_Results/LUSC_EA/broad_significance_results.txt', sep='\t')

# df of alteration freq with annotation
df_sq_Amp=data.frame(Arm=sq_AA$Arm, AA_Freq=sq_AA$Amp.frequency, EA_Freq=sq_EA$Amp.frequency, AA_qvalue=sq_AA$Amp.q.value, EA_qvalue=sq_EA$Amp.q.value,
                        CNV_Type='Amp')
df_sq_Del=data.frame(Arm=sq_AA$Arm, AA_Freq=sq_AA$Del.frequency, EA_Freq=sq_EA$Del.frequency, AA_qvalue=sq_AA$Del.q.value, EA_qvalue=sq_EA$Del.q.value,
                        CNV_Type='Del')
df_adeno_Amp=data.frame(Arm=adeno_AA$Arm, AA_Freq=adeno_AA$Amp.frequency, EA_Freq=adeno_EA$Amp.frequency, AA_qvalue=adeno_AA$Amp.q.value, EA_qvalue=adeno_EA$Amp.q.value,
                        CNV_Type='Amp')
df_adeno_Del=data.frame(Arm=adeno_AA$Arm, AA_Freq=adeno_AA$Del.frequency, EA_Freq=adeno_EA$Del.frequency, AA_qvalue=adeno_AA$Del.q.value, EA_qvalue=adeno_EA$Del.q.value,
                        CNV_Type='Del')
# Adding recurrence infomation for Adeno
df_adeno_Amp$Sig = 'Recurrent in None'
df_adeno_Amp$Sig[df_adeno_Amp$AA_qvalue<0.05 & df_adeno_Amp$EA_qvalue<0.05] = 'Recurrent in Both'
df_adeno_Amp$Sig[df_adeno_Amp$AA_qvalue<0.05 & df_adeno_Amp$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_adeno_Amp$Sig[df_adeno_Amp$AA_qvalue>0.05 & df_adeno_Amp$EA_qvalue<0.05] = 'Recurrent in EA Only'
df_adeno_Del$Sig = 'Recurrent in None'
df_adeno_Del$Sig[df_adeno_Del$AA_qvalue<0.05 & df_adeno_Del$EA_qvalue<0.05] = 'Recurrent in Both'
df_adeno_Del$Sig[df_adeno_Del$AA_qvalue<0.05 & df_adeno_Del$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_adeno_Del$Sig[df_adeno_Del$AA_qvalue>0.05 & df_adeno_Del$EA_qvalue<0.05] = 'Recurrent in EA Only'

df_adeno=rbind(df_adeno_Amp, df_adeno_Del)
# Adding recurrence infomation for sq
df_sq_Amp$Sig = 'Recurrent in None'
df_sq_Amp$Sig[df_sq_Amp$AA_qvalue<0.05 & df_sq_Amp$EA_qvalue<0.05] = 'Recurrent in Both'
df_sq_Amp$Sig[df_sq_Amp$AA_qvalue<0.05 & df_sq_Amp$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_sq_Amp$Sig[df_sq_Amp$AA_qvalue>0.05 & df_sq_Amp$EA_qvalue<0.05] = 'Recurrent in EA Only'
df_sq_Del$Sig = 'Recurrent in None'
df_sq_Del$Sig[df_sq_Del$AA_qvalue<0.05 & df_sq_Del$EA_qvalue<0.05] = 'Recurrent in Both'
df_sq_Del$Sig[df_sq_Del$AA_qvalue<0.05 & df_sq_Del$EA_qvalue>0.05] = 'Recurrent in AA Only'
df_sq_Del$Sig[df_sq_Del$AA_qvalue>0.05 & df_sq_Del$EA_qvalue<0.05] = 'Recurrent in EA Only'

df_sq=rbind(df_sq_Amp, df_sq_Del)

# Formatting the dfs
df=rbind(cbind(df_adeno, Hist='LUAD'), cbind(df_sq, Hist='LUSC'))
colnames(df)[7]='Significance'
df_NCIMD=df
########################
# Figure 2B #in TCGA
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
# Vector form of Figure 2
########################
setEPS()
postscript('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Figure2_vecForm.eps',
     width = 24, height = 13.3)
plot(AX)
dev.off()
