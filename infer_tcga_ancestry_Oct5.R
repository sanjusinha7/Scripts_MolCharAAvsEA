#####################################################
###MYCUSTOM functions
#####################################################
###Following are custom functions I've made to avoid repetitively writing same things again - after three years of my PhD :p 
myhead<-function(x){
  x[1:min(5, nrow(x)), 1:min(5, ncol(x))]
}
err_handle<-function(x){ tryCatch(x, error=function(e){NA}) }
test_topbottomQuantile<-function(var1=UCD_score_M1, var2_tobetiled=unlist(Exp_metab[1,]), which_tail='g', numofQuantiles=3){
  length(var1); length(var2_tobetiled)
  var2_tobetiled=xtile(var2_tobetiled, numofQuantiles)
  ##Only keep top/bottom quantile
  var1= var1[var2_tobetiled== min(var2_tobetiled) | var2_tobetiled== max(var2_tobetiled)]
  var2_tobetiled= var2_tobetiled[var2_tobetiled== min(var2_tobetiled) | var2_tobetiled== max(var2_tobetiled)] 
  c(sig=err_handle(wilcox.test(var1 ~ factor(var2_tobetiled), alternative=which_tail)$p.value),
    eff_size=err_handle(diff(aggregate(var1, by=list(var2_tobetiled), mean)[,2]) ) )
}
hypergeometric_test_for_twolists<-function(test_list, base_list, global, lowertail=FALSE) {
  #If lowertail=FALSE - we calculate the probability for enrichment
  length(base_list)
  adj_base_list=global[na.omit(match(base_list, global))]
  Matched_list=test_list[!is.na(match(test_list, adj_base_list))]
  phyper(length(Matched_list)-1, length(adj_base_list), length(global)- length(adj_base_list), length(test_list), lower.tail=lowertail)
}
fdrcorr<-function(test_list){p.adjust(test_list, method = 'fdr')}

# dftemp=data.frame(prob$samples, prob$types, prob$stage)
# sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUSC',], function(x) length(x))[1:3,3])/sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUSC',], function(x) length(x))[,3])
# sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUAD',], function(x) length(x))[1:4,3])/sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUAD',], function(x) length(x))[,3])
# split(aggregate(mat$X ~ mat$stage_trimmed+mat$hist, data=dftemp, function(x) length(x)))
# 

#####################################################
###Goal: Using PCA- infer ancestry
#####################################################
# library(GenomicDataCommons)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# require(devtools)
# install_github("Bioconductor/GenomicDataCommons")
require(GenomicDataCommons)
library(magrittr)
library(ggfortify)
library(cluster)
require('lfda')
require('e1071')
require('caret')
require('pROC')
require("factoextra")
require(FactoMineR)
require(parallel)
require(data.table)
require(tictoc)
load('/data/sinhas8/ISLE-master/data/TCGA.RData')
prob$race=factor(prob$race)
levels(prob$race)[c(3,6)]=c('AA', 'EA')
TCGAtranslateID = function(file_ids, legacy = TRUE) {
  info = files(legacy = legacy) %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  # The mess of code below is to extract TCGA barcodes
  # id_list will contain a list (one item for each file_id)
  # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
  id_list = lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  # so we can later expand to a data.frame of the right size
  barcodes_per_file = sapply(id_list,length)
  # And build the data.frame
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}
read_genotype<-function(infunc_folderName=nested_folderList[1], cancer_type, flag){
  setwd(paste('/data/sinhas8/tcga_genotype/retrieval_', cancer_type, '_201909/', infunc_folderName, sep=''))
  print(flag)
  getwd()
  ##Search the genptype file by the name of the algorithm in the filename used to generate
  fileName=list.files()[grep('birdseed',list.files())]
  SNP=fread(fileName, sep='\t')
  SNP[,2]
}
unique(prob$types)
####Function to infer ancestry
infer_ancestry<-function(cancer_type='KICH', test_run=F, samples_used_in_test=20){
  #####################################################
  ###Initialize Step 0: Choosing the cancer type
  #####################################################
  print(paste('Initialized for', cancer_type))
  print(paste('started for',cancer_type))
  folderName=paste('/data/sinhas8/tcga_genotype/retrieval_',cancer_type,'_201909/', sep='')
  setwd(folderName)
  nested_folderList=list.files()
  tcgaid_genotype=TCGAtranslateID(nested_folderList)
  ###Only consider matched to prob samples
  # nested_folderList_matched=nested_folderList[!is.na(match(substring(tcgaid_genotype[,2], 0, 12), prob$samples))]
  # tcgaid_genotype_matched=tcgaid_genotype[!is.na(match(substring(tcgaid_genotype[,2], 0, 12), prob$samples)),]
  tcgaid_genotype=data.frame(tcgaid_genotype)
  tcgaid_genotype$sample=substring(tcgaid_genotype$submitter_id, 0, 12)
  tcgaid_genotype$type=substring(tcgaid_genotype[,2], 14, )
  tcgaid_genotype$whether_normal=as.numeric(substring(tcgaid_genotype[,2], 14, 15))>=10
  tcgaid_genotype=tcgaid_genotype[tcgaid_genotype$whether_normal,]
  tcgaid_genotype=tcgaid_genotype[match(unique(tcgaid_genotype$sample), tcgaid_genotype$sample),]
  tcgaid_genotype=tcgaid_genotype[!is.na(match(tcgaid_genotype$sample, prob$samples[prob$types==cancer_type])),]
  #####################################################
  ###started Step 0-demo init: Create Additional TCGA dataset
  #####################################################
  print(paste('started Step 0-demo init'))
  tcgaid_genotype$race=prob$race[match(tcgaid_genotype$sample, prob$samples)]
  tcgaid_genotype$race=factor(tcgaid_genotype$race)
  tcgaid_genotype$stage=prob$stage[match(tcgaid_genotype$sample, prob$samples)]
  tcgaid_genotype$sex=prob$sex[match(tcgaid_genotype$sample, prob$samples)]
  tcgaid_genotype$age=prob$age[match(tcgaid_genotype$sample, prob$samples)]
  nested_folderList=as.character(tcgaid_genotype$file_id)
  print(length(nested_folderList))
  print(table(tcgaid_genotype$race))
  #####################################################
  ###started Step 1-infer genotype: Create Genotype
  #####################################################
  print(paste('started Step 1-infer genotype'))
  ##Create genotype matrix
  if(test_run){
    tic()
    genotyping_list=mclapply(1:samples_used_in_test, function(x) err_handle(read_genotype(nested_folderList[x],cancer_type, x)), mc.cores = detectCores())
    genotyping_mat=do.call(cbind, genotyping_list)
    tcgaid_genotype=tcgaid_genotype[1:samples_used_in_test,]
    toc()
  } else{
    tic()
    genotyping_list=mclapply(1:length(nested_folderList), function(x) err_handle(read_genotype(nested_folderList[x], cancer_type, x)))
    genotyping_mat=do.call(cbind, genotyping_list)
    toc()
  }
  
  colnames(genotyping_mat)=as.character(tcgaid_genotype[,3])
  genotyping_mat=genotyping_mat[-1,]
  saveRDS(genotyping_mat, paste('/data/sinhas8/tcga_ancestry_inference/genotyping_mat_complete', cancer_type, '.RDS', sep=''))
  genotyping_mat_df=as.data.frame(genotyping_mat)
  genotyping_mat_df_var=apply(genotyping_mat_df, 1, var)
  Variance_Threshold=sum(tcgaid_genotype$race=='AA', na.rm = T)/sum(tcgaid_genotype$race=='EA'|tcgaid_genotype$race=='AA', na.rm = T)
  genotyping_mat_df=genotyping_mat_df[which(genotyping_mat_df_var>Variance_Threshold),]
  AA_pop=which(tcgaid_genotype$race=='AA'); EA_pop=which(tcgaid_genotype$race=='EA')
  AA_pop_genoDist=sapply(1:nrow(genotyping_mat_df), function(x) table(factor(unlist(genotyping_mat_df[x,AA_pop]), levels=c('0','1','2')) ))
  EA_pop_genoDist=sapply(1:nrow(genotyping_mat_df), function(x) table(factor(unlist(genotyping_mat_df[x,EA_pop]), levels=c('0','1','2')) ))
  table_pop=t(rbind(AA=AA_pop_genoDist, EA=EA_pop_genoDist))
  saveRDS(table_pop, paste('/data/sinhas8/tcga_ancestry_inference/table_pop', cancer_type, '.RDS', sep=''))
  #####################################################
  ##started Step 2-pca: Choose Race associated SNVs
  #####################################################
  print(paste('started Step 2-pca'))
  K=0.50
  AA_assoc_fullallele=apply(table_pop, 1, function(x) x[3]> sum(x[1:3])*K )
  EA_assoc_fullallele=apply(table_pop, 2, function(x) x[6]> sum(x[4:6])*K )
  dim(genotyping_mat_df[which(AA_assoc_fullallele | EA_assoc_fullallele),])
  mat=genotyping_mat_df[which((AA_assoc_fullallele | EA_assoc_fullallele)),]
  mat=apply(mat, 1, as.numeric)
  saveRDS(mat, paste('/data/sinhas8/tcga_ancestry_inference/geno_mat_', cancer_type, '.RDS', sep=''))
  forced_pca_try1=prcomp(mat, rank=2)
  saveRDS(forced_pca_try1, paste('/data/sinhas8/tcga_ancestry_inference/forced_pca_try1', cancer_type, '.RDS', sep=''))
  #####################################################
  ###started Step 3-plotting: Plotting
  #####################################################
  print(paste('started Step 3-plotting'))
  tiff(paste('/data/sinhas8/tcga_ancestry_inference/pca_plot', cancer_type, '.RDS', sep=''))
  autoplot(forced_pca_try1, data = tcgaid_genotype, colour = 'race', frame=F)+theme_classic(base_size = 20)
  dev.off()
  #####################################################
  ###started Step 4-infer ancestry via svm: Classification using SVM
  #####################################################
  print(paste('started Step 4-infer ancestry via svm'))
  forsvm_raw=data.frame(PC=forced_pca_try1$x, race=tcgaid_genotype$race)
  # forsvm=forsvm_raw[is.na(forsvm_raw$race) | forsvm_raw$race=='AA' | forsvm_raw$race=='EA',]
  forsvm=forsvm_raw[forsvm_raw$race=='AA' | forsvm_raw$race=='EA',]
  forsvm$race=factor(as.character(forsvm$race))
  forsvm=forsvm[,c(1,3)]
  svmfit=svm(race~., data=forsvm, kernel='linear', scale=F)
  predictions=factor(as.numeric(svmfit$decision.values<0), labels = c('AA','EA'))
  length(predictions)
  results=confusionMatrix(predictions,
                          na.omit(forsvm$race))
  results
  tcgaid_genotype=tcgaid_genotype[which(tcgaid_genotype$race=='AA' | tcgaid_genotype$race=='EA'),]
  tcgaid_genotype$inferred_ancestry=predictions
  saveRDS(tcgaid_genotype, paste('/data/sinhas8/tcga_ancestry_inference/results_p1/',cancer_type, '_infAncestry.RDS', sep=''))
  tcgaid_genotype
}

setwd('/data/sinhas8/tcga_genotype')
cancer_type_list=list.dirs(recursive = F)
cancer_type_list=sapply(cancer_type_list, function(x) strsplit(x, '_')[[1]][2])
cancer_type_list=cancer_type_list[match(levels(factor(prob$types)), cancer_type_list)]
length(cancer_type_list)
inferred_ancestry=sapply(cancer_type_list[-c(3, 11,28)], function(x) err_handle(infer_ancestry(cancer_type=x, test_run=F, samples_used_in_test = 10)) )
saveRDS(inferred_ancestry, '/data/sinhas8/tcga_ancestry_inference/results_p1/Complete_Inferred_Ancestry.RDS')

