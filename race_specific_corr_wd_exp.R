#Function to test the significnace of correlation btw a CNV change of a gene adn expression
#########################
##Step 0: Define functions
#########################
race_specific_corr_wd_exp<-function(GeneList=adeno_list, Cancer_Type='adeno', demo=mat){
	#Files Needed
	#CNV matrices of a given cancer type
	CNV_AA=read.csv(paste('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/AA_',Cancer_Type, '/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
	CNV_EA=read.csv(paste('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/EA_',Cancer_Type, '/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]

	##Preprocessing::
	#removing  the '[' or ']' signs from Gene-names,
	GeneList=sapply(GeneList, function(x) gsub('\\]','',gsub('\\[','',x)))
	#Final CNV matrix
	CNV=cbind(CNV_AA, CNV_EA)
	#removing duplicated rows
	CNV=CNV[-which(duplicated(sapply(rownames(CNV), function(x) strsplit(x, '\\|')[[1]][1]))),]	
	rownames(CNV)=sapply(rownames(CNV), function(x) strsplit(x, '\\|')[[1]][1]) 
	colnames(CNV)= sapply(colnames(CNV), function(x) gsub('X','',x))

	##Matching rows
	common_genes_2all= Reduce(intersect, list(GeneList, rownames(Exp), rownames(CNV))) 
	demo$names_modified=gsub('-','.',demo$names)

	#Making colnames uniform
	#remove "X" from colanmes due to names inconsistency.
	colnames(CNV)= demo$acc[na.omit(match(colnames(CNV), demo$names_modified))]
	colnames(Exp)= sapply(colnames(Exp), function(x) gsub('T','',strsplit(x, '_')[[1]][1]))

	Matched_Cols = Reduce(intersect,list(colnames(Exp), colnames(CNV)))
	Exp_matched  = Exp[common_genes_2all, Matched_Cols]
	CNV_matched  = CNV[common_genes_2all, Matched_Cols]

	#calculating wilcoxon rank sum test
	wil_mat=sapply(1:nrow(Exp_matched), function(i) tryCatch(wilcox.test(unlist(as.numeric(Exp_matched[i,])) ~ as.factor(xtile(unlist(CNV_matched[i,]),2)), alternative='l')$p.value, error=function(e){NA}))

	names(wil_mat) = common_genes_2all
	wil_mat[order(wil_mat)]
}

