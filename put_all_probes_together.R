##Put all probes together
#Phase1
library(cluster)
require('ggfortify')
require('factoextra')
require(tidyr)
require(parallel)
require("ggbiplot")
Phase1='/Users/sinhas8/Project_Chromotrypsis/2.Data/Oncosan_Phase1/tsv/'
Phase2='/Users/sinhas8/Project_Chromotrypsis/2.Data/New_OncoScan/tsv/'
annotation=

final_mat<-function(foldername){
  setwd(foldername)
  filelist=list.files()
  genotype_filelist=filelist[grep('Calls',filelist)]
  l<-mclapply(genotype_filelist , function(file) { 
    r <- cbind(file, read.csv(file, skip = 16, header = T, sep='\t')[,c(2,3,5)])
    r
  }, mc.cpres=detectCores())
  df= do.call(rbind , l)
  mat=spread(df[,1:3],ProbeSetName, Call)
  forced_mat=spread(df[,c(1,2,4)],ProbeSetName, ForcedCall)
  list(mat, forced_mat)
}

test1=final_mat(Phase1)
names(test1)=c('mat','forced_mat')
test2=final_mat(Phase2)
names(test2)=c('mat','forced_mat')
mat=rbind(test1$mat, test2$mat)
forced_mat=rbind(test1$forced_mat, test2$forced_mat)
dim(mat); dim(forced_mat)
###Pre-processing to fit standards
##6 = AA genotype##7 = BB genotype##8 = AB genotype##11 = NoCall
saveRDS(mat,'/Users/sinhas8/Project_Chromotrypsis/2.Data/genotype_matrix.RDS')
saveRDS(forced_mat,'/Users/sinhas8/Project_Chromotrypsis/2.Data/forced_genotype_matrix.RDS')

##
mat=readRDS('/Users/sinhas8/Project_Chromotrypsis/2.Data/genotype_matrix.RDS')
forced_mat=readRDS('/Users/sinhas8/Project_Chromotrypsis/2.Data/forced_genotype_matrix.RDS')
rownames(forced_mat)=forced_mat$file

###PCA_analysis
#mat[mat==6]=0;mat[mat==7]=2;mat[mat==8]=1;mat[mat==11]=NA
forced_mat[forced_mat==6]=0
forced_mat[forced_mat==7]=2
forced_mat[forced_mat==8]=1
forced_mat[forced_mat==11]=NA

#####matching
demo=read.csv('/users/sinhas8/Project_Chromotrypsis/2.Data/Corrected_NCIMD_HRD_by_LOH_and_GI.csv')
demo=demo[demo$hist=='adeno' | demo$hist =='sq',]
demo$hist=as.character(demo$hist)
demo$hist=factor(demo$hist)
levels(demo$hist)=c('LUAD', 'LUSC')

rownames(forced_mat)[-grep('_',rownames(forced_mat))]=sapply(gsub('duplicate','',rownames(forced_mat)[-grep('_',rownames(forced_mat))]), function(x) strsplit(x, '\\.')[[1]][1])
rownames(forced_mat)[grep('Recentering', rownames(forced_mat))]='1_581T,Genotyping,Calls.tsv'
rownames(forced_mat)[grep('_',rownames(forced_mat))]=sapply(rownames(forced_mat)[grep('_',rownames(forced_mat))],
                                                            function(x) strsplit(strsplit(x, ',')[[1]][1], '_')[[1]][2] )
rownames(forced_mat)=gsub('[Tt]','',rownames(forced_mat))
forced_mat=forced_mat[!is.na(match(rownames(forced_mat), demo$acc)),]

forced_mat=forced_mat[,-1]
forced_mat[1:5, 1:5]
forced_pca_try1=prcomp(forced_mat, rank=1)
clustering=clara(forced_pca_try1$x, k=2)
#forced_pca_try2=prcomp(forced_mat, center = T, scale. = T, rank=3)

autoplot(forced_pca_try1, data = demo, colour = 'race', frame=TRUE)


