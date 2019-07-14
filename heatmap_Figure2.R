###
t1=df1_pp[grepl('Del',df1_pp[,2]) | grepl('Del',df1_pp[,3]),c(1,6,7)]
colnames(t1)=c('GeneName','Freq_EA','Freq_AA')
t2=df1_pp[!(grepl('Del',df1_pp[,2]) | grepl('Del',df1_pp[,3])),c(1,4,5)]
colnames(t2)=c('GeneName','Freq_EA','Freq_AA')
df2plot=rbind(t1, t2)
rownames(df2plot)=df2plot[,1]
df2plot=df2plot[,-1]

require('NMF')
tiff('./Project_Chromotrypsis/prep_final_figures/heatmap_1_scaled.tif', 400, 400)
aheatmap(df2plot, Rowv = F, Colv = NA, 
         fontsize = 10, width = 300, height = 300,
         annRow = data.frame(Event_Type=c(rep('Del', nrow(df2plot)-2), 
                                                        rep('Amp', 2))),
         scale='row')
dev.off()
