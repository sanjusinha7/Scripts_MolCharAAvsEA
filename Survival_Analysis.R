##Survival Analysis
##Required libraries
require(survival)

#List of genes associated with Expression in Squamous carcinoma and Adenocarcinoma
AssocwdExp_Genes_adeno=readRDS('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/Assocwd_Exp_Adeno.RDS')
AssocwdExp_Genes_sq=readRDS('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/Assocwd_Exp_Squamous.RDS')


name_it<-function(Cancer_Type, GeneList=AssocwdExp_Genes_adeno){
	CNV_AA=read.csv(paste('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/AA_',Cancer_Type, '/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
	CNV_AA=CNV_AA[-which(duplicated(sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
	rownames(CNV_AA)=sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]) 

	CNV_EA=read.csv(paste('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/EA_',Cancer_Type, '/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
	CNV_EA=CNV_EA[-which(duplicated(sapply(rownames(CNV_EA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
	rownames(CNV_EA)=sapply(rownames(CNV_EA), function(x) strsplit(x, '\\|')[[1]][1]) 
	
	CNV=cbind(CNV_AA, CNV_EA)
	#Demographics
	mat=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/demo_and_clones.csv')
	mat=mat[grep(Cancer_Type,mat$hist),]

	GeneList=names(which(GeneList<0.1))
	CNV_FF=CNV[match(GeneList, rownames(CNV)),]

	colnames(CNV_FF)=sapply(colnames(CNV_FF), function(x) substring(x, 2))
	mat_FF=mat[na.omit(match(sapply(colnames(CNV_FF), function(x) gsub('-','\\.',x)),  sapply(mat$names, function(x) gsub('-','\\.', x)))),]
	CNV_SF=CNV_FF[,colnames(CNV_FF)[!is.na(match(sapply(colnames(CNV_FF), function(x) gsub('-','\\.',x)),  sapply(mat, function(x) gsub('-','\\.', x))))]]

	qq=apply(CNV_SF, 1, function(x) summary(coxph(Surv(mat_FF$survival, mat_FF$lungcancer_death_all_years) ~ x + mat_FF$stage_trimmed+ mat_FF$race))$coefficients[1,c(1,5)])

	previous_df=read.csv(paste('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Freq_GeneLevel_for_reccurentin_', Cancer_Type, '.csv', sep=''))
	rr=data.frame(previous_df[match(colnames(qq), previous_df$X),], Survival_cox_=t(qq))
	
	rr[rr$Survival_cox_.Pr...z..<0.1,]

}

adeno_cox_results=name_it(Cancer_Type='adeno', GeneList=AssocwdExp_Genes_adeno)
sq_cox_results=name_it(Cancer_Type='sq', GeneList=AssocwdExp_Genes_sq)

write.csv(adeno_cox_results, '/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/intermediate_adeno.csv', quote=F, row.names=F)
write.csv(sq_cox_results, '/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/intermediate_sq.csv', quote=F, row.names=F)

