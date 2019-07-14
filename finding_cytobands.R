##Finding Cytobands given a gene and a cancer type.
AssocwdExp_Genes_adeno=as.data.frame(AssocwdExp_Genes_adeno)
AssocwdExp_Genes_sq=as.data.frame(AssocwdExp_Genes_sq)

finding_cytobands<-function(AssocwdExp_Genes=AssocwdExp_Genes_adeno, AA=adeno_AA, EA=adeno_EA){
	AA_ID=sapply(rownames(AssocwdExp_Genes), function(x) as.character(AA$Descriptor)[grep(paste(x,'\\|',sep=''), as.character(AA$Genes))])
	EA_ID=sapply(rownames(AssocwdExp_Genes), function(x) as.character(EA$Descriptor)[grep(paste(x,'\\|',sep=''), as.character(EA$Genes))])

	AA_ID_percentage=sapply(rownames(AssocwdExp_Genes), function(x) as.character(AA$Frequency)[grep(paste(x,'\\|',sep=''), as.character(AA$Genes))])
	EA_ID_percentage=sapply(rownames(AssocwdExp_Genes), function(x) as.character(EA$Frequency)[grep(paste(x,'\\|',sep=''), as.character(EA$Genes))])

	AA_ID_percentage=AA_ID_percentage[sapply(AA_ID_percentage,length)>0]
	AA_ID_percentage=sapply(AA_ID_percentage, function(x) paste(x, collapse='\\|'))	
	EA_ID_percentage=EA_ID_percentage[sapply(EA_ID_percentage,length)>0]
	EA_ID_percentage=sapply(EA_ID_percentage, function(x) paste(x, collapse=','))	

	AA_ID=AA_ID[sapply(AA_ID,length)>0]
	AA_ID=sapply(AA_ID, function(x) paste(x, collapse='\\|'))	
	EA_ID=EA_ID[sapply(EA_ID,length)>0]
	EA_ID=sapply(EA_ID, function(x) paste(x, collapse=','))	

	AssocwdExp_Genes$AA_cytoband = ''
	AssocwdExp_Genes$EA_cytoband = ''

	AssocwdExp_Genes$AA_percentage = 0
	AssocwdExp_Genes$EA_percentage = 0

	##
	AssocwdExp_Genes$AA_cytoband[!is.na(match(rownames(AssocwdExp_Genes), names(AA_ID[sapply(AA_ID,length)>0]) ))]=AA_ID
	AssocwdExp_Genes$EA_cytoband[!is.na(match(rownames(AssocwdExp_Genes), names(EA_ID[sapply(EA_ID,length)>0]) ))]=EA_ID
	AssocwdExp_Genes$AA_percentage[!is.na(match(rownames(AssocwdExp_Genes), names(AA_ID_percentage) ))]=AA_ID_percentage
	AssocwdExp_Genes$EA_percentage[!is.na(match(rownames(AssocwdExp_Genes), names(EA_ID_percentage) ))]=EA_ID_percentage

	AssocwdExp_Genes
}

adeno= finding_cytobands(AssocwdExp_Genes=AssocwdExp_Genes_adeno, AA=adeno_AA, EA=adeno_EA)
sq= finding_cytobands(AssocwdExp_Genes=AssocwdExp_Genes_sq, AA=sq_AA, EA=sq_EA)

##
write.csv(adeno, '/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Corr_wd_Expression/Association_test_adeno.csv', quote=F)
write.csv(sq, '/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Corr_wd_Expression/Association_test_Sq.csv', quote=F)

