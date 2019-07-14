##
sq_Ass=readRDS('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Corr_wd_Expression/Association_test_sq.RDS')
adeno_Ass=readRDS('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Corr_wd_Expression/Association_test_adeno.RDS')

from_PLOS_SQ=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/plos2016paper_squamous.csv')
from_PLOS_ad=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/plos2016paper_adeno.csv')

##Testing ENrichemnt
testing_enrichment<-function(from_PLOS, from_GISTIC_corr_Exp, Sig_Threshold=0.1, Cancer_Type='adeno'){
	##
	base_list = as.character(from_PLOS[which(as.numeric(as.character(from_PLOS$FDR))<0.01),1])
#	universe  = rownames(from_GISTIC_corr_Exp)

	CNV_AA=read.csv(paste('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/AA_',Cancer_Type, '/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
	CNV_EA=read.csv(paste('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/EA_',Cancer_Type, '/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]

	#Final CNV matrix
	CNV=cbind(CNV_AA, CNV_EA)
	CNV=CNV[-which(duplicated(sapply(rownames(CNV), function(x) strsplit(x, '\\|')[[1]][1]))),]	
	rownames(CNV)=sapply(rownames(CNV), function(x) strsplit(x, '\\|')[[1]][1]) 

	universe=rownames(CNV)
	test_list = rownames(from_GISTIC_corr_Exp)[which(from_GISTIC_corr_Exp[,1] < Sig_Threshold)]	
	print(data.frame(length(base_list), length(universe),length(test_list)))
	hypergeometric_test_for_twolists(test_list, base_list, universe)

}

##Hypergeometric test
hypergeometric_test_for_twolists<-function(test_list, base_list, global, lowertail=TRUE) {
	adj_base_list=global[na.omit(match(base_list, global))]
	Matched_list=test_list[!is.na(match(test_list, adj_base_list))]
	phyper(length(test_list)-length(Matched_list), length(global)- length(adj_base_list), length(adj_base_list), length(test_list), lower.tail=lowertail)
#	Matched_list
}

##Calling the above enrichement test function
#Robustness: Test the below function iwht various threshold of Sig_Thresh = c(0.01, 0.05, 0.1, 0.2)
testing_enrichment(from_PLOS=from_PLOS_SQ, from_GISTIC_corr_Exp=sq_Ass, Sig_Threshold=0.1, Cancer_Type='sq')
testing_enrichment(from_PLOS=from_PLOS_ad, from_GISTIC_corr_Exp=adeno_Ass, Sig_Threshold=0.1, Cancer_Type='adeno')

##############################################################################################################
#####################################Integrating vaious columns###############################################
##############################################################################################################
sq_Ass=readRDS('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Corr_wd_Expression/Association_test_sq.RDS')
adeno_Ass=readRDS('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Corr_wd_Expression/Association_test_adeno.RDS')


df_integration<-function(sq_Ass, Cancer_Type){
	sq_Ass_Sig=sq_Ass[which(sq_Ass[,1]<0.1),]

	##CNV dataset
	CNV_AA=read.csv(paste('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/AA_',Cancer_Type, '/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
#	CNV=cbind(CNV_AA, CNV_EA)
	CNV_AA=CNV_AA[-which(duplicated(sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
	rownames(CNV_AA)=sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]) 

	CNV_EA=read.csv(paste('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/EA_',Cancer_Type, '/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
	#Final CNV matrix
	CNV_EA=CNV_EA[-which(duplicated(sapply(rownames(CNV_EA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
	rownames(CNV_EA)=sapply(rownames(CNV_EA), function(x) strsplit(x, '\\|')[[1]][1]) 

	GeneList=rownames(sq_Ass_Sig)	
	df_AMP= data.frame(Amp_EA_freq=apply(CNV_EA[GeneList,]>0.1,1,sum)/ncol(CNV_EA), Amp_AA_freq=apply(CNV_AA[GeneList,]>0.1,1,sum)/ncol(CNV_AA))
	df_DEL= data.frame(Del_EA_freq=apply(CNV_EA[GeneList,]< -0.1,1,sum)/ncol(CNV_EA), Del_AA_freq=apply(CNV_AA[GeneList,]< -0.1,1,sum)/ncol(CNV_AA))
	df=cbind(df_AMP, df_DEL)



}
