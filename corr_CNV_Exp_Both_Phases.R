##Script to get gene level copy number from segment level:: (TCGA level 3 from TCGA level 2 file.)
##libraries
#########################
##Step 0: Libraries
#########################
require(data.table)
require(statar)
#########################
##Step 0: Files needed
#########################
setwd('/Users/sinhas8/Project_Chromotrypsis/')
source('3.Tools/race_specific_corr_wd_exp.R')
Exp=read.csv('/Users/sinhas8/Downloads/RSEM_CPM_TMM_counts.txt', sep='\t')	#Exp matrix
#Demographics
mat=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/demo_and_clones.csv')
##Below are frequency files.
LUAD_AA=read.csv('4.Results/gistic_TCGA/LUAD_AA/freq_of_cytobands_wd_Genes.csv')
LUAD_EA=read.csv('4.Results/gistic_TCGA/LUAD_EA/freq_of_cytobands_wd_Genes.csv')
LUSC_AA=read.csv('4.Results/gistic_TCGA/LUSC_AA/freq_of_cytobands_wd_Genes.csv')
LUSC_EA=read.csv('4.Results/gistic_TCGA/LUSC_EA/freq_of_cytobands_wd_Genes.csv')

#########################
##Pre-Processing
#########################
Exp=as.matrix(Exp)
Exp_genelist= sapply(Exp[,1], function(x) strsplit(x, '\\|')[[1]][2])
Exp=Exp[,-1]
rownames(Exp)=Exp_genelist
Exp=Exp[,grep('T',colnames(Exp))]

LUAD_AA_list=  unique(unlist(sapply(LUAD_AA$Genes, function(x) unlist(strsplit(as.character(x), '\\| ')))))
LUAD_EA_list=  unique(unlist(sapply(LUAD_EA$Genes, function(x) unlist(strsplit(as.character(x), '\\| ')))))
LUSC_AA_list   =  unique(unlist(sapply(LUSC_AA$Genes, function(x) unlist(strsplit(as.character(x), '\\| ')))))
LUSC_EA_list   =  unique(unlist(sapply(LUSC_EA$Genes, function(x) unlist(strsplit(as.character(x), '\\| ')))))

##Taking a union of the two population recurrent changes
LUAD_geneList=union(LUAD_AA_list, LUAD_EA_list)
LUSC_geneList=   union(LUSC_AA_list, LUSC_EA_list)

##
#########################
##Step 0: Define functions
#########################
Step1<-function(Cancer_Type='LUAD'){
  if(Cancer_Type=='LUSC'){
    AssocwdExp_Genes=NA
    GeneList=LUSC_geneList; EA=LUSC_EA;AA=LUSC_AA
  } else{
    AssocwdExp_Genes=NA
    GeneList=LUAD_geneList; EA=LUAD_EA;AA=LUAD_AA
    }
  GeneList=sapply(GeneList, function(x) gsub('\\[|\\]','',x))
  ##
  EA_Genes=lapply(as.character(EA$Genes), function(x) strsplit(gsub('\\[|\\]|\\s','',x),'\\|'))
  AA_Genes=lapply(as.character(AA$Genes), function(x) strsplit(gsub('\\[|\\]|\\s','',x),'\\|'))
  
  Aberrance_Type_EA=sapply(GeneList, function(x) paste(substring(as.character(EA$Unique.Name[which(unlist(sapply(EA_Genes, function(y) length(na.omit(match(y[[1]], x)))>0 )))]), 1, 3), collapse='|') )
  Aberrance_Type_AA=sapply(GeneList, function(x) paste(substring(as.character(AA$Unique.Name[which(unlist(sapply(AA_Genes, function(y) length(na.omit(match(y[[1]], x)))>0 )))]), 1, 3), collapse='|') )
  
  cytoband_EA = sapply(GeneList, function(x) paste(as.character(EA$Descriptor[which(unlist(sapply(EA_Genes, function(y) length(na.omit(match(y[[1]], x)))>0 )))]), collapse='|') )
  cytoband_AA = sapply(GeneList, function(x) paste(as.character(AA$Descriptor[which(unlist(sapply(AA_Genes, function(y) length(na.omit(match(y[[1]], x)))>0 )))]), collapse='|') )
  
  ##Significance
  Sig_EA = sapply(GeneList, function(x) paste(as.character(EA$q.values[which(unlist(sapply(EA_Genes, function(y) length(na.omit(match(y[[1]], x)))>0 )))]), collapse='|') )
  Sig_AA = sapply(GeneList, function(x) paste(as.character(AA$q.values[which(unlist(sapply(AA_Genes, function(y) length(na.omit(match(y[[1]], x)))>0 )))]), collapse='|') )
  
  ##CNV dataset for gene 
  CNV_AA=read.csv(paste('4.Results/gistic_TCGA/',Cancer_Type, '_AA','/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
  CNV_AA=CNV_AA[-which(duplicated(sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
  rownames(CNV_AA)=sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]) 
  
  CNV_EA=read.csv(paste('4.Results/gistic_TCGA/',Cancer_Type, '_EA','/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
  CNV_EA=CNV_EA[-which(duplicated(sapply(rownames(CNV_EA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
  rownames(CNV_EA)=sapply(rownames(CNV_EA), function(x) strsplit(x, '\\|')[[1]][1]) 
  
  ##
  df_AMP = data.frame(Amp_EA_freq=apply(CNV_EA[GeneList,]> 0.5,1,sum)/ncol(CNV_EA), Amp_AA_freq=apply(CNV_AA[GeneList,]> 0.5,1,sum)/ncol(CNV_AA))
  df_DEL = data.frame(Del_EA_freq=apply(CNV_EA[GeneList,]< -0.5,1,sum)/ncol(CNV_EA), Del_AA_freq=apply(CNV_AA[GeneList,]< -0.5,1,sum)/ncol(CNV_AA))
  
  #df1=data.frame(Signigicance=AssocwdExp_Genes[GeneList], Aberrance_AA=Aberrance_Type_adeno_AA, Aberrance_EA=Aberrance_Type_adeno_EA, cyto_AA=cytoband_adeno_AA, cyto_EA=cytoband_adeno_EA, df_AMP, df_DEL)
  
  df1=data.frame(Sig_EA, Sig_AA,  Aberrance_EA=Aberrance_Type_EA, Aberrance_AA=Aberrance_Type_AA, cyto_EA=cytoband_EA, cyto_AA=cytoband_AA, df_AMP, df_DEL)
  df1$Exp_Assoc_Sig=NA
  
  df1$Exp_Assoc_Sig[na.omit(match(names(AssocwdExp_Genes), rownames(df1)))] =AssocwdExp_Genes[!is.na(match(names(AssocwdExp_Genes), rownames(df1)))]
  ##
  df1
  
}
genes_aff_survival<-function(genelist_EA_Amp, race='EA', Type){
  CNV=read.csv(paste('/home/sinhas8/Both_',race,'/all_data_by_genes.txt', sep=''), sep='\t', row.names=1)
  CNV=CNV[,-c(1:2)]
  CNV_genelist=sapply(rownames(CNV), function(x) strsplit(x, '\\|')[[1]][1]) 
  colnames(CNV)=sapply(colnames(CNV), function(x) gsub('X','',x))
  CNV=CNV[na.omit(match(genelist_EA_Amp, CNV_genelist)),]
  colnames(CNV)=mat$acc[na.omit(match(colnames(CNV), mat$names_modified))]
  
  temp_mat=mat[na.omit(match(colnames(CNV), mat$acc)),]	
  
  cox_reg=do.call(rbind,lapply(1:nrow(CNV), function(x) summary(coxph(Surv(temp_mat$survival, temp_mat$lungcancer_death_all_years)~ unlist(CNV[x,])))[[7]]))
  rownames(cox_reg)=genelist_EA_Amp
  cox_reg=cox_reg[,c(1,5)]
  cox_reg=cox_reg[order(cox_reg[,2]),]
  write.csv(cox_reg ,paste('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/', 'Survival',Type,race,'.csv', sep=''))
}
#Enrichment test
hypergeometric_test_for_twolists<-function(test_list, base_list, global, lowertail=TRUE) {
  adj_base_list=global[na.omit(match(base_list, global))]
  Matched_list=test_list[!is.na(match(test_list, adj_base_list))]
  phyper(length(test_list)-length(Matched_list), length(global)- length(adj_base_list), length(adj_base_list), length(test_list), lower.tail=lowertail)
}
#########################
##Step 0: Callfunctions for TCGA
#########################
df_LUSC=Step1('LUSC')
df_LUAD=Step1('LUAD')

write.csv(df_LUSC, '4.Results/gistic_TCGA/Freq_LUSC', quote=F)
write.csv(df_LUAD, '4.Results/gistic_TCGA/Freq_LUAD', quote=F)
#########################
##Step 2: Survival PArt
#########################
df1_Amp=df1[grepl('Amp', df1$Aberrance_AA) | grepl('Amp', df1$Aberrance_EA),]
df1_Del=df1[grepl('Del', df1$Aberrance_AA) | grepl('Del', df1$Aberrance_EA),]

genelist_Amp_EA=names(which(p.adjust(unlist(CorrwdExp_Genes_Amp_EA[,1]), method='fdr')<0.1))
genelist_Del_EA=names(which(p.adjust(unlist(CorrwdExp_Genes_Del_EA[,1]), method='fdr')<0.1))
genelist_Amp_AA=names(which(p.adjust(unlist(CorrwdExp_Genes_Amp_AA[,1]), method='fdr')<0.1))
genelist_Del_AA=names(which(p.adjust(unlist(CorrwdExp_Genes_Del_AA[,1]), method='fdr')<0.1))
  