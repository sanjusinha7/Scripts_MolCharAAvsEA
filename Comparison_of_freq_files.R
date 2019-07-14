##Function to get the events having atleast one known cancer genes.
COSMIC=read.csv('/home/sinhas8/Downloads/COSMIC.tsv', sep='\t')
events_wd_COSMICgenes<- function(Filename=''){
	adeno_AA=read.csv(Filename)	
	OncoGenes=COSMIC$Gene.Symbol[grep('oncogene',COSMIC$Role.in.Cancer)]
	TSG=COSMIC$Gene.Symbol[grep('TSG',COSMIC$Role.in.Cancer)]

	COSMIC_Onco_Mapping_regions=sapply(OncoGenes, function(x) grep(paste(x,'\\|',sep=''), adeno_AA$Genes[grep('Amp', adeno_AA$Unique.Name)]))

	COSMIC_TSG_Mapping_regions=sapply(TSG, function(x) grep(paste(x,'\\|',sep=''), adeno_AA$Genes[grep('Del', adeno_AA$Unique.Name)]))

	Total_Del=length(adeno_AA$Genes[grep('Del', adeno_AA$Unique.Name)])
	Total_Amp=length(adeno_AA$Genes[grep('Amp', adeno_AA$Unique.Name)])
	Known_OncoGenes_regions=length(unique(unlist(COSMIC_Onco_Mapping_regions)))
	Known_TSG_regions=length(unique(unlist(COSMIC_TSG_Mapping_regions)))

	data.frame(Total_Del, Known_TSG_regions, Total_Amp, Known_OncoGenes_regions)
#	COSMIC$Gene.Symbol[which(sapply(COSMIC_Gene_Mapping, function(x) length(x)>0))]
}

adeno_AA=events_wd_COSMICgenes('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/AA_adeno/freq_of_cytobands_wd_Genes_adeno_AA.csv')
adeno_EA=events_wd_COSMICgenes('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/EA_adeno/freq_of_cytobands_wd_Genes_EA_adeno.csv')
sq_AA=events_wd_COSMICgenes('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/AA_sq/freq_of_cytobands_wd_Genes_AA_sq.csv')
sq_EA=events_wd_COSMICgenes('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/EA_sq/freq_of_cytobands_wd_Genes_EA_sq.csv')


###
comparison_of_known_onco<- function(Filename='/home/sinhas8/GISTIC_Results_byraceandhistology/AA_adeno/freq_of_cytobands_wd_Genes_adeno_AA.csv'){
	adeno_AA=read.csv(Filename)
	COSMIC_Gene_Mapping=sapply(COSMIC$Gene.Symbol, function(x) grep(paste(x,'\\|',sep=''), adeno_AA$Genes))
	data.frame(length(unique(unlist(COSMIC_Gene_Mapping))), nrow(adeno_AA))
#	COSMIC$Gene.Symbol[which(sapply(COSMIC_Gene_Mapping, function(x) length(x)>0))]

	adeno_AA[which(sapply(COSMIC_Gene_Mapping, function(x) length(x)>0)),]
}


##Venn Diagram
for_Venn_Diagram<-function(f1, f2){
	file1=read.csv(f1)
	file2=read.csv(f2)
	#For Deletion first
	g1=data.frame(AA=file1$Descriptor[grep('Del',file1$Unique.Name)], EA=file2$Descriptor[grep('Del',file2$Unique.Name)])
	g2=cbind(AA=file1[grep('Amp',file1$Unique.Name),]$Descriptor, EA=file2[grep('Amp',file2$Unique.Name),]$Descriptor)


	intersection_Del=length(intersect(file1_Del$Descriptor, file2_Del$Descriptor))
	area1_Del=sum(is.na(match(file1_Del$Descriptor, file2_Del$Descriptor)))
	area2_Del=sum(is.na(match(file2_Del$Descriptor, file1_Del$Descriptor)))
	
	intersection_Amp=length(intersect(file1_Amp$Descriptor, file2_Amp$Descriptor))
	area1_Amp=sum(is.na(match(file1_Amp$Descriptor, file2_Amp$Descriptor)))
	area2_Amp=sum(is.na(match(file2_Amp$Descriptor, file1_Amp$Descriptor)))
	list(data.frame(intersection_Del, area1_Del, area2_Del), data.frame(intersection_Amp, area1_Amp, area2_Amp))

	
}


adeno_venn=for_Venn_Diagram(f1='/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/AA_adeno/freq_of_cytobands_wd_Genes_adeno_AA.csv', f2='/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/EA_adeno/freq_of_cytobands_wd_Genes_EA_adeno.csv')

sq_venn=for_Venn_Diagram('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/AA_sq/freq_of_cytobands_wd_Genes_AA_sq.csv', f2='/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/EA_sq/freq_of_cytobands_wd_Genes_EA_sq.csv')

