##DEV PHASE
#THEME:: The Alchemist :: The secret of life, though, is to fall seven times and to get up eight times.”

#Loading reaquired files
CNV=readRDS('/cbcb/project2-scratch/sanju/Chromotrypsis/2.Data/CG0314_GSA_Initial_standard_cnv_report.RDS')
CNV=CNV[order(CNV$'Sample ID', CNV$Chr, as.numeric(CNV$Start)),]
CNV_List= split(CNV, CNV$'Sample ID')

############################################################################################################################################################################################
#Hypothesis 1:: Signal-background Hypothesis::The alteration events on a Chromotripsis Chromosome(CC) should be higher than the background events. 
############################################################################################################################################################################################
source('/cbcb/project2-scratch/sanju/Chromotrypsis/3.Tools/Events_Enrichment_Test.R')
Test_H1=lapply(CNV_List, function(x)H1_Events_Enrichment(x))
############################################################################################################################################################################################
#Hypothesis 2 test:: Single events wwould be enriched in the chromosomes which has been through chromotripsis.:: phyper(q, m, n, k)
############################################################################################################################################################################################
source('/cbcb/project2-scratch/sanju/Chromotrypsis/3.Tools/Oscillatory_Events_Enrichement_Test.R')
###!!!!!!!!!!Change the THreshold to 0.5 in final analysis to check the results change.!!!!!!!!
Test_H2= lapply(CNV_List, Oscillatory_Enrichment_Test)

############################################################################################################################################################################################
#Hypothesis 3:: Breakpoint CLustering:: Breakpoints should be closer than background distance.
############################################################################################################################################################################################
source('/cbcb/project2-scratch/sanju/Chromotrypsis/3.Tools/BPC_Test_Illumina.R')
Test_H3 = lapply(CNV_List, BPC_Test)


#Hyperparamters
Hypergeom_Threshold=0.5
CHTP_Events_Count_Threshold=10
P.Value_Threshold=0.01
With_H1=0

##Hypothesis 1 and 2 together. 
Result_Hyperparameters<-function(Hypergeom_Threshold, CHTP_Events_Count_Threshold, P.Value_Threshold, With_H1=1){
	if(With_H1==1){
		Enrichment_Inference=sapply(H12, function(x) rownames(x)[x[,1] & x[,2]<Hypergeom_Threshold])
	}
	else{
		Enrichment_Inference=sapply(H12, function(x) rownames(x)[x[,2]<Hypergeom_Threshold])
	}
	H3_Inference=sapply(Test_H3, function(x) names(x)[sapply(x, function(y) y['P.Value']<P.Value_Threshold & y['CHTP_Events_Count']>CHTP_Events_Count_Threshold)] )
	Result_Try1=lapply(1:length(CNV_List), function(x) na.omit(H3_Inference[[x]][match(Enrichment_Inference[[x]], H3_Inference[[x]])]))
	list(length(which(sapply(Result_Try1, length)>0)), sort(table(unlist(Enrichment_Inference))), sort(table(na.omit(unlist(Result_Try1)))))
}

Result_Hyperparameters(0.5, 2, 0.1, 0)


#It's the possibility of having a dream come true that makes life interesting.
