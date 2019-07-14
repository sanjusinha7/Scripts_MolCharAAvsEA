#This script is comprised of functions to define a feature of Chromotripsis called:: Breakpoint Cluster.
source('/cbcb/project2-scratch/sanju/Chromotrypsis/3.Tools/Events_Enrichment_Test.R')
source('/cbcb/project2-scratch/sanju/Chromotrypsis/3.Tools/BPC_Test.R')
source('/cbcb/project2-scratch/sanju/Chromotrypsis/3.Tools/Oscillatory_Events_Enrichement_Test.R')


#We need to define a few functions are varibale to use this script.

############################################################################################################################################################################################
#Analyzing distribution of events!!
############################################################################################################################################################################################
for (i in Samples_wd_Chtp){
	#Order:: q_list, m, n, k_list
	print(c(as.character(CNV_Type2_List[[i]][1,2]), i))
	q_list=sapply(CHTP_Distribution[[i]], function(x) sum(x, na.rm=T))
	m=sum(unlist(CHTP_Distribution[[i]]), na.rm=T)
	n=sum(unlist(Distribution[[i]]), na.rm=T)-m
	k_list=sapply(Distribution[[i]], sum)

	print(sapply(CHTP_Distribution[[i]], function(x) sum(x, na.rm=T)))
	print(sum(unlist(CHTP_Distribution[[i]]), na.rm=T))
	print(sum(unlist(Distribution[[i]]), na.rm=T)-m)
	print(sapply(Distribution[[i]], sum))
#	print(sum(p.adjust(mapply(function(q, k) phyper(q-1, m, n, k, lower.tail=F), q_list, k_list)[1:22], method='fdr')<0.2))
}

