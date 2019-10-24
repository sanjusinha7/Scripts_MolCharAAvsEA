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
require(parallel)
require(data.table)
require(tictoc)
read_genotype<-function(infunc_folderName=nested_folderList[1]){
  paste('/Users/sinhas8/Project_Chromotrypsis/2.Data/retrieval_UCEC_201909/', infunc_folderName, sep='')
  ##Search the genptype file by the name of the algorithm in the filename used to generate
  fileName=list.files()[grep('birdseed',list.files())]
  SNP=fread(fileName, sep='\t')
  SNP[,2]
}
folderName='/Users/sinhas8/Project_Chromotrypsis/2.Data/retrieval_UCEC_201909'
setwd(folderName)
nested_folderList=list.files()
length(nested_folderList)

tic()
UCEC_genotyping_list=lapply(1:100, function(x) err_handle(read_genotype(nested_folderList[x])))
toc()
object.size(UCEC_genotyping_list)/1024/1024/1024
UCEC_genotyping=do.call(cbind, UCEC_genotyping_list)
myhead(UCEC_genotyping)
#####################################################
###Mapping to samples
#####################################################
lcd /Volumes/CDSL/TCGA/
put /Volumes/CDSL/TCGA/retrieval_BR*
