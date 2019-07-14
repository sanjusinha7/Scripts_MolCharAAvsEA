##Hierachial Clustering o AA and AA

adeno_AA=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/AA_adeno/broad_values_by_arm.txt', sep='\t', row.names=1)
adeno_EA=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/EA_adeno/broad_values_by_arm.txt', sep='\t', row.names=1)

sq_AA=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/AA_sq/broad_values_by_arm.txt', sep='\t', row.names=1)
sq_EA=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/EA_sq/broad_values_by_arm.txt', sep='\t', row.names=1)

#Annotation=c(rep('AA',ncol(AA)), rep('EA', ncol(EA)))


hcluster <- function(Event_Matrix, Type='Complete'){
	Event_Matrix_dist=dist(Event_Matrix)
	clusters=hclust(Event_Matrix_dist)
	tiff(paste('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/', Type, 'Broad_Level_Events_Corr.tiff', sep=''))
	plot(clusters)
	dev.off()
}


hcluster(adeno_AA, 'adeno_AA')
hcluster(adeno_EA, 'adeno_EA')
hcluster(sq_AA, 'sq_AA')
hcluster(sq_EA, 'sq_EA')


