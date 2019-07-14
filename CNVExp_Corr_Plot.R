##libraries
require(data.table)
require(statar)

#########################
##Step 0: Pre-Processing
#########################
##functions needed
getwd()
setwd('/Users/sinhas8/Project_Chromotrypsis/')
source('./3.Tools/race_specific_corr_wd_exp.R')
##Expression
Exp=cbind(read.csv('2.Data/Exp_AA_AD.csv'),
             read.csv('2.Data/Exp_EA_AD.csv'),
             read.csv('2.Data/Exp_AA_SC.csv'),
             read.csv('2.Data/Exp_EA_SC.csv'))
rownames(Exp)=make.names(Exp[,1], unique = T)
Exp=Exp[,-1]
Exp=Exp[,-grep('X',colnames(Exp))]
colnames(Exp)=sapply(colnames(Exp), function(x) gsub("[^0-9]", "", x))
##CNV dataset for gene 
setwd('/Users/sinhas8/Project_Chromotrypsis/')
Cancer_Type='LUSC'
CNV_AA=read.csv(paste('./4.Results/GISTIC_results/Trimmed_GISTIC_Results/',Cancer_Type, '_AA','/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
CNV_AA=CNV_AA[-which(duplicated(sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
rownames(CNV_AA)=sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]) 
CNV_EA=read.csv(paste('./4.Results/GISTIC_results/Trimmed_GISTIC_Results/',Cancer_Type, '_EA','/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
CNV_EA=CNV_EA[-which(duplicated(sapply(rownames(CNV_EA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
rownames(CNV_EA)=sapply(rownames(CNV_EA), function(x) strsplit(x, '\\|')[[1]][1]) 
CNV1=cbind(CNV_AA, CNV_EA)
Cancer_Type='LUAD'
CNV_AA=read.csv(paste('./4.Results/GISTIC_results/Trimmed_GISTIC_Results/',Cancer_Type, '_AA','/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
CNV_AA=CNV_AA[-which(duplicated(sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
rownames(CNV_AA)=sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]) 
CNV_EA=read.csv(paste('./4.Results/GISTIC_results/Trimmed_GISTIC_Results/',Cancer_Type, '_EA','/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
CNV_EA=CNV_EA[-which(duplicated(sapply(rownames(CNV_EA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
rownames(CNV_EA)=sapply(rownames(CNV_EA), function(x) strsplit(x, '\\|')[[1]][1]) 
CNV2=cbind(CNV_AA, CNV_EA)
CNV=cbind(CNV1, CNV2)
colnames(CNV)[grep('_',colnames(CNV))]= 
  sapply(gsub('X','',colnames(CNV)[grep('_',colnames(CNV))]), function(x) strsplit(x,'_')[[1]][2])
colnames(CNV)=gsub('X','',colnames(CNV))
colnames(CNV)[1]=substring(colnames(CNV)[1], 1, nchar(colnames(CNV)[1])-10)
colnames(CNV)=sapply(sapply(colnames(CNV), function(x) strsplit(x, '\\.')[[1]][1]),
                     function(x) gsub("[^0-9]", "", x))
#########################
##Step 1: Matched Matrix
#########################
CNV_matched=CNV[,na.omit(match(colnames(Exp), colnames(CNV)))]
Exp_matched=Exp[,!is.na(match(colnames(Exp), colnames(CNV)))]
CNV_matched=CNV_matched[na.omit(match(rownames(Exp), rownames(CNV))),]
Exp_matched=Exp_matched[!is.na(match(rownames(Exp), rownames(CNV))),]
Exp_matched_scaled=apply(Exp_matched, 1, scale)
CNV_matched_scaled=apply(CNV_matched, 1, scale)
#########################
##Step 2: Matched Matrix
#########################
dim(CNV_matched)
corr_CNV=sapply(seq(nrow(CNV_matched)), function(x) 
  unlist(cor.test(unlist(CNV_matched[x,]), unlist(Exp_matched[x,]))[c(3, 4)]))
corr_CNV[1:2, 1:10]

require(ggExtra)
df2plot=data.frame(Exp_Z=c(Exp_matched_scaled), 
                   CNV_Z=c(CNV_matched_scaled),
                   CNV  =c(as.matrix(CNV_matched)),
                   Event_Type='') 
df2plot$Event_Type=as.character(df2plot$Event_Type)
df2plot$Event_Type[df2plot$CNV> 2.9]='Amp'
df2plot$Event_Type[df2plot$CNV< 1.1]='Del'
cor.test(df2plot$Exp_Z, df2plot$CNV_Z)
head(df2plot)
Panel2<-ggplot(data=df2plot, aes(x=Exp_Z, y=CNV_Z))+
  geom_point(data=df2plot, aes(x=Exp_Z, y=CNV_Z, color=Event_Type), alpha=0.01)+
  stat_smooth(method = 'lm')+
  annotate("text", label = "Pearson Rho=0.21, p<1.6E-16", x = 2, y = 6, size=8)+
  theme_bw(base_size = 20)+
  labs(x=' Expression z-score', y=' CNV z-score')+
  scale_color_manual(values=c("black",'Red','Green'))+
  theme()
plot(p)
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/corr_plot.tif')
ggMarginal(p, type = "density", size = 4)
dev.off()

