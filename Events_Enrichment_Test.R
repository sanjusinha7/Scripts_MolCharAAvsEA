#This script is comprised of functions to define a feature of Chromotripsis called:: Events Enrichment.

############################################################################################################################################################################################
#Hypothesis 1:: Signal-background Hypothesis::The alteration events on a Chromotripsis Chromosome(CC) should be higher than the background events. 
############################################################################################################################################################################################
H1_Events_Enrichment<-function(CNV_Test){
	Temp=split(CNV_Test, CNV_Test$Chromosome)
	Threshold=mean(sapply(Temp, function(x)sum(table(x$CN.State))))
	as.matrix(sapply(sapply(Temp, function(x) sum(table(x$CN.State))), function(y) y>Threshold))
}

