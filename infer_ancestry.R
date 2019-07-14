library(ggfortify)
library(cluster)
require('lfda')
require('e1071')
require('caret')
require('pROC')
install.packages(c("FactoMineR", "factoextra"))
forced_mat=readRDS('/Users/sinhas8/Project_Chromotrypsis/2.Data/forced_genotype_matrix.RDS')
mat=read.csv('/users/sinhas8/Project_Chromotrypsis/2.Data/Corrected_NCIMD_HRD_by_LOH_and_GI.csv')
mat=mat[mat$hist=='adeno' | mat$hist =='sq',]
mat$hist=as.character(mat$hist)
mat$hist=factor(mat$hist)
levels(mat$hist)=c('LUAD', 'LUSC')

#
rownames(forced_mat)=forced_mat$file
forced_mat=forced_mat[,-1]

rownames(forced_mat)=gsub('duplicate','',rownames(forced_mat))
rownames(forced_mat)[-grep('_',rownames(forced_mat))]=sapply(rownames(forced_mat)[-grep('_',rownames(forced_mat))], function(x) 
  strsplit(x, '\\.')[[1]][1])
rownames(forced_mat)[grep('_',rownames(forced_mat))][1]=substring(rownames(forced_mat)[grep('_',rownames(forced_mat))][1], 1, 4)
rownames(forced_mat)[grep('_',rownames(forced_mat))]=sapply(rownames(forced_mat)[grep('_',rownames(forced_mat))], 
                                                            function(x) strsplit(strsplit(x, ',')[[1]][1], '_')[[1]][2]  )
rownames(forced_mat)=gsub('[Tt]','',rownames(forced_mat))
forced_mat=forced_mat[!is.na(match(rownames(forced_mat), mat$acc)),]
dim(forced_mat)

forced_mat[1:5, 1:5]
forced_pca_try1=prcomp(forced_mat, rank=2)

tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/PCA_inferenceRace.tiff')
autoplot(forced_pca_try1, data = mat, colour = 'race')+theme_classic(base_size = 20)
dev.off()

autoplot(clara(forced_pca_try1$x,k=2), data=mat, color='race')
#forced_pca_try2=prcomp(forced_mat, center = T, scale. = T, rank=3)

pam_results=lfda(as.matrix(forced_pca_try1$x[,1]),r=2)
table(pam_results$clustering, mat$race)

###Cluster using SVM
head(forsvm)
forsvm=data.frame(PC=forced_pca_try1$x, race=mat$race)
forsvm=forsvm[,c(1,3)]
tail(forsvm)
svmfit=svm(race~., data=forsvm, kernel='linear', scale=F)
plot(svmfit, data=forsvm)
predictions=factor(as.numeric(svmfit$decision.values<0), labels = c('AA','EA'))
mat$inferred_ancestry=predictions
results=confusionMatrix(predictions,
                mat$race)
results
