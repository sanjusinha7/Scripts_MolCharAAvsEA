##Testign Signature 3 across race
err_handle<-function(x){ tryCatch(x, error=function(e){NA}) }
TCGAsig=read.csv('/Users/sinhas8/Downloads/signature_profile_sample.txt', sep='\t')
head(TCGAsig)
levels(TCGAsig$race)[c(3,6)]=c('AA','EA')
TCGAsig=TCGAsig[TCGAsig$race=='AA' | TCGAsig$race=='EA',]
cancer_type_list=levels(TCGAsig$Project_Name)
signature_mapping= as.numeric(gsub('Signature.','',levels(TCGAsig$Signature)))

Diff_across_race<-function(cancer_type=levels(TCGAsig$Project_Name)[43],
                           Sig_Num=3, 
                           pancan=FALSE){
  which_signature=levels(TCGAsig$Signature)[order(signature_mapping)[Sig_Num]]
  if(pancan){
    mat=TCGAsig[which(
      grepl('TCGA',TCGAsig$Project_Name) &  TCGAsig$Signature == which_signature
    ),]
  } else{
  mat=TCGAsig[which(
    TCGAsig$Project_Name == cancer_type & TCGAsig$Signature == which_signature
    ),]
  }
  wilcox.test(mat$Contribution ~ factor(as.character(mat$race)), alternative='g')$p.value
}

TCGA_projects=levels(TCGAsig$Project_Name)[grep('TCGA',levels(TCGAsig$Project_Name))]

TCGA_Signature13_cancerTYpe=data.frame(TCGA_projects, P=sapply(1:73, function(x)
  err_handle(Diff_across_race(cancer_type=levels(TCGAsig$Project_Name)[x],
                              Sig_Num=13, pancan=FALSE)))[grep('TCGA',
                                                              levels(TCGAsig$Project_Name))])
write.csv(TCGA_Signature12_cancerTYpe, '/Users/sinhas8/TCGA_Signature12_cancerTYpe.csv')

TCGA_Signature3_cancerTYpe=na.omit(TCGA_Signature3_cancerTYpe)
sum(TCGA_Signature3_cancerTYpe$P<0.5)/nrow(TCGA_Signature3_cancerTYpe)

SigFor_HRD=sapply(1:30, function(x)
  err_handle(Diff_across_race(cancer_type=levels(TCGAsig$Project_Name)[1], Sig_Num=x, pancan=TRUE)))

Pancan_Sigdf=data.frame(Sig=levels(TCGAsig$Signature)[order(signature_mapping)[1:30]],SigFor_HRD)
Pancan_Sigdf=Pancan_Sigdf[order(Pancan_Sigdf$SigFor_HRD),]
Pancan_Sigdf$FDR_SigFor_HRD=p.adjust(Pancan_Sigdf$SigFor_HRD, method='fdr')

table(TCGAsig$race[TCGAsig$Project_Name==levels(TCGAsig$Project_Name)[43]])/30
################
##Signature 3
################
dim(TCGAsig)
TCGAsig=TCGAsig[grep('TCGA',TCGAsig$Project_Name),]
sig3_TCGA=TCGAsig[TCGAsig$Signature=="Signature.3",]



