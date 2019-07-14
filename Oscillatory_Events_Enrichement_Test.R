#This script is comprised of functions to define a feature of Chromotripsis called:: Oscillatory Events enrichment.

############################################################################################################################################################################################
#Hypothesis 2 test:: Single events wwould be enriched in the chromosomes which has been through chromotripsis.:: phyper(q, m, n, k)
############################################################################################################################################################################################

#Events distribution
Events_Distribution <- function(CNV_Test){
	Temp_H1=split(CNV_Test, CNV_Test$Chromosome)
	lapply(Temp_H1, function(x) table(x$CN.State))
}

#RDS:: Results Determing step :: Oscillation events distribution for hypergeometric test
CHTP_Events_Distribution <- function(CNV_Test){
	Temp_H1=split(CNV_Test, CNV_Test$Chromosome)
	lapply(Temp_H1, function(x) table(x$CN.State)[c('1', '3')])		#Condsidering 2.5 and 1.5 as 1 and 3 to not minimize any true negative.
}


Oscillatory_Enrichment_Test<-function(Input){
	Distribution=Events_Distribution(Input)
	CHTP_Distribution=CHTP_Events_Distribution(Input)
	q_list=sapply(CHTP_Distribution, function(x) sum(x, na.rm=T))
	m=sum(unlist(CHTP_Distribution), na.rm=T)
	n=sum(unlist(Distribution), na.rm=T)-m
	k_list=sapply(Distribution, sum)


	mapply(function(q, k) phyper(q-1, m, n, k, lower.tail=F), q_list, k_list)
}

