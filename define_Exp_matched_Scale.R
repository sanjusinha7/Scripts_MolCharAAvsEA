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
