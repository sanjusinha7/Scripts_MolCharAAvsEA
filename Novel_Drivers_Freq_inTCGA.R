#########
##Testing for Driver frequency patterns in TCGA
#########
prob=readRDS('../sinhas8/Downloads/TCGA_withMut.RDS')
prob$race=factor(prob$race, labels = c('AI','AS','AA','HL','NH','EA') )
GOI=as.character(Pot_Novel_Drivers$Gene.Name)
Freq_across_Population<-function(GeneName, cancer_type='LUSC'){
  geneID=which(prob$genes==GeneName)
  sample_OI=which(prob$types==cancer_type)
  df=as.data.frame(aggregate(prob$scna[geneID,sample_OI], list(prob$race[sample_OI]) , 
                             function(x) c(Del=sum(x< -0.1)/length(x), Amp=sum(x> 0.1)/length(x)) ))
  df=data.frame(Del=df[,2][,1], Amp=df[,2][,2])
  c(df[2,1], df[4,1], df[2,2], df[4,2])
}
df_GOI=data.frame(t(apply(Pot_Novel_Drivers, 1, function(x) 
  Freq_across_Population(as.character(x[1]),
                         cancer_type=as.character(x[13])) )))
rownames(df_GOI)=GOI
colnames(df_GOI)=c('AA_Del', 'EA_Del', 'AA_Amp', 'EA_Amp')
df_GOI$hist='LUAD'
#df_GOI$Type=c(rep('Amp',3), 'Del')
df_GOI$Type='Del'
df=df_GOI
df_Amp=df[df$Type=='Amp',c(3,4, 5:6),]
colnames(df_Amp)[1:2]=c('AA_Freq','EA_Freq')
df_Del=df[df$Type=='Del',c(1,2, 5:6),]
colnames(df_Del)[1:2]=c('AA_Freq','EA_Freq')
df=rbind(df_Amp, df_Del)
df$GeneName=rownames(df)

Figure6A=ggplotGrob(
  ggplot(df, aes(x=AA_Freq, y=EA_Freq, label = GeneName))+
    geom_point(position=position_jitter(h=0.01,w=0.01), size=2.5)+ 
    geom_abline(intercept = 0, slope = 1, color='blue')+
    geom_text_repel(size=8, show.legend = FALSE)+ 
    facet_grid(hist~Type)+
    theme_bw(base_size = 24)+
    labs(x='Frequency in AA', y='Frequency in EA')+
    xlim(range=c(0,0.75))+
    ylim(range=c(0,0.75))
  # scale_color_manual(values=c("blue", "orange", "red", "black"))+ 
  # guides(color=guide_legend(nrow=2, byrow=TRUE, override.aes = list(size = 3) ))+
  # #   guides(color=guide_legend()+
  # theme(legend.position="top", legend.text=element_text(size=20),
  #       legend.title=element_text(size=22))
)
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Novel_Drivers_TCGA_corrected.tif',
     height = 800, width = 800)
plot(Figure6A)
dev.off()
