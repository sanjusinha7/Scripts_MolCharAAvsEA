incosmic<-function(df=adeno_AA){
	df$cosmic_matched_Genes=''
	df_list=sapply(df$Genes, function(x) strsplit(as.character(x),'\\| ')[[1]])

	Mapping = unlist(sapply(as.character(cosmic$Gene.Symbol), function(x) 
                       which(sapply(df_list, function(y) sum(!is.na(match(x, y)))>0 ))))
	Map_list = lapply(split(data.frame(Mapping), Mapping), rownames)
	df$cosmic_matched_Genes[as.numeric(names(Map_list))] = 
           sapply(Map_list, function(x) paste(x, collapse='| '))

	df
}

adeno_AA_incosmic=incosmic(adeno_AA)
adeno_EA_incosmic=incosmic(adeno_EA)
sq_AA_incosmic=incosmic(sq_AA)
sq_EA_incosmic=incosmic(sq_EA)

write.csv(adeno_AA_incosmic,'/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/adeno_AA_incosmic.csv')
write.csv(adeno_EA_incosmic,'/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/adeno_EA_incosmic.csv')
write.csv(sq_AA_incosmic,'/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/sq_AA_incosmic.csv')
write.csv(sq_EA_incosmic,'/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/sq_EA_incosmic.csv')


#######
unique(sq_AA_incosmic$Descriptor[which(sq_AA_incosmic$cosmic_matched_Genes!='')])


