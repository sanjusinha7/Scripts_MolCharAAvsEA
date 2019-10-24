#########
##Testing for Driver frequency patterns in TCGA
#########
prob=readRDS('/Users/sinhas8/Downloads/TCGA_withMut.RDS')
prob$race=factor(prob$race, labels = c('AI','AS','AA','HL','NH','EA') )

Driver_GeneList=unique(Drivers$GeneName)
Driver_GeneList=c(Driver_GeneList, 'PTEN')
Freq_across_Population<-function(GeneName, cancer_type='LUSC'){
  geneID=which(prob$genes==GeneName)
  sample_OI=which(prob$types==cancer_type)
  df=as.data.frame(aggregate(prob$scna[geneID,sample_OI], list(prob$race[sample_OI]) , 
               function(x) c(Del=sum(x< -0.1)/length(x), Amp=sum(x> 0.1)/length(x)) ))
  df=data.frame(Del=df[,2][,1], Amp=df[,2][,2])
  c(df[2,1], df[4,1], df[2,2], df[4,2])
}
df_LUSC=data.frame(t(sapply(as.character(Driver_GeneList), function(x) 
  Freq_across_Population(x, cancer_type='LUSC') )))
rownames(df_LUSC)=Driver_GeneList
colnames(df_LUSC)=c('AA_Del', 'EA_Del', 'AA_Amp', 'EA_Amp')
df_LUSC$Type=c(rep('Amp',5), rep('Del',nrow(df_LUSC)-5))
df_LUSC$GeneName=rownames(df_LUSC)
df_LUSC$hist='LUSC'

df_LUAD=data.frame(t(sapply(Driver_GeneList, 
                            function(x) Freq_across_Population(x, cancer_type='LUAD') )))
rownames(df_LUAD)=Driver_GeneList
colnames(df_LUAD)=c('AA_Del', 'EA_Del', 'AA_Amp', 'EA_Amp')
df_LUAD$Type=c(rep('Amp',5), rep('Del',nrow(df_LUAD)-5))
df_LUAD$GeneName=rownames(df_LUAD)
df_LUAD$hist='LUAD'

df=rbind(df_LUSC, df_LUAD)
df_Amp=df[df$Type=='Amp',c(3,4, 5:7),]
colnames(df_Amp)[1:2]=c('AA_Freq','EA_Freq')
df_Del=df[df$Type=='Del',c(1,2, 5:7),]
colnames(df_Del)[1:2]=c('AA_Freq','EA_Freq')
df=rbind(df_Amp, df_Del)
Figure6A=ggplotGrob(
  ggplot(df, aes(x=AA_Freq, y=EA_Freq, label = GeneName))+
    geom_point(position=position_jitter(h=0.01,w=0.01), size=2.5)+ 
    geom_abline(intercept = 0, slope = 1, color='blue')+
    geom_text_repel(size=8, show.legend = FALSE)+ 
    facet_grid(hist~Type)+
    theme_bw(base_size = 24)+
    labs(x='Frequency in AA', y='Frequency in EA')+
    ggtitle('TCGA Lung cancer Drivers Alteration Freq across race')
    # scale_color_manual(values=c("blue", "orange", "red", "black"))+ 
    # guides(color=guide_legend(nrow=2, byrow=TRUE, override.aes = list(size = 3) ))+
    # #   guides(color=guide_legend()+
    # theme(legend.position="top", legend.text=element_text(size=20),
    #       legend.title=element_text(size=22))
)
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Figure6A_TCGA_May8th.tif',
     height = 1200, width = 1200)
plot(Figure6A)
dev.off()
