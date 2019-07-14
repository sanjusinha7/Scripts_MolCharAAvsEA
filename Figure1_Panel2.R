## Plot of frequency distribution of recurretly aberrant genes.
##*Follow up of corr_CNV_Exp_Both_Phases.R*

df=df1_Del[,c('Del_AA_freq','Del_EA_freq')]
ScatterPlot<-ggplot(df, aes(x=Del_AA_freq, y=Del_EA_freq ))+
	geom_point(position=position_jitter(h=0.01,w=0.01))+
	geom_abline(intercept = 0, slope = 1, color='blue')+
	lims( x=c(0,max(df$Del_AA_freq, df$Del_EA_freq)), y=c(0,max(df$Del_AA_freq, df$Del_EA_freq)) )

xdensity <- ggplot(df, aes(Del_AA_freq)) + 
  geom_density(alpha=.5)

ydensity <- ggplot(df, aes(Del_EA_freq)) + 
  geom_density(alpha=.5)
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(), 
   panel.border = element_blank(),
   panel.background = element_blank(),
   axis.title.x = element_blank(),
   axis.title.y = element_blank(),
   axis.text.x = element_blank(), 
   axis.text.y = element_blank(),
   axis.ticks = element_blank()
     )

tiff('/home/sinhas8/Projects/Project_Chromotrypsis/Results_New/adeno_Del_4th.tiff')
grid.arrange(xdensity, blankPlot, ScatterPlot, ydensity, 
        ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4) )
dev.off()

