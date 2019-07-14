##DEV PHASE
#THEME:: The Alchemist :: The secret of life, though, is to fall seven times and to get up eight times.”

#Loading Files Required
FilenameQC=read.csv('/cbcb/project2-scratch/sanju/Chromotrypsis/2.Data/Filename_QC.txt',sep='\t')
CNV1=read.csv('/cbcb/project2-scratch/sanju/Chromotrypsis/2.Data/CN_State_of_Segments_CHAS_wdlocation.txt', sep='\t')


############################################################################################################################################################################################
#Hypothesis 1:: Signal-background Hypothesis::The alteration events on a Chromotripsis Chromosome(CC) should be higher than the background events. 
############################################################################################################################################################################################
source('/cbcb/project2-scratch/sanju/Chromotrypsis/3.Tools/Events_Enrichment_Test.R')

#Only-Diploid_Samples
Diploid_Samples=(FilenameQC[,22]==2)
Diploid_Samples[is.na(Diploid_Samples)]=FALSE
Diploid_Tumour_Samples=(FilenameQC$File[Diploid_Samples])[-c(grep('n\\.',FilenameQC$File[Diploid_Samples]), grep('N\\.',FilenameQC$File[Diploid_Samples]))]

#Ordering the sample based on Chr location
CNV1=CNV1[order(CNV1$Full.Location),]

#Breakdown into a list
CNV_Type2=CNV1[!is.na(match(CNV1$File, Diploid_Tumour_Samples)),]
CNV_Type2_List=split(CNV_Type2, CNV_Type2$File)
CNV_Type2_List= CNV_Type2_List[sapply(CNV_Type2_List, function(x) dim(x)[1]>0)]

#CNV_Type2_Dip_List=CNV_Type2_List[match(FilenameQC$File[which(FilenameQC$tuscanPloidy..CHP.Summary.<3)], names(CNV_Type2_List))]

H1_Events_Enrichment=lapply(CNV_Type2_List, function(x)H1_Events_Enrichment(x))

############################################################################################################################################################################################
#Hypothesis 2 test:: Single events wwould be enriched in the chromosomes which has been through chromotripsis.:: phyper(q, m, n, k)
############################################################################################################################################################################################
source('/cbcb/project2-scratch/sanju/Chromotrypsis/3.Tools/Oscillatory_Events_Enrichement_Test.R')

Distribution=lapply(CNV_Type2_List, function(x) Events_Distribution(x) )
CHTP_Distribution=lapply(CNV_Type2_List, function(x) CHTP_Events_Distribution(x) )


H2_Results=c();Id_Flag=1
for (i in 1:length(CNV_Type2_List)) {
	q_list=sapply(CHTP_Distribution[[i]], function(x) sum(x, na.rm=T))
	m=sum(unlist(CHTP_Distribution[[i]]), na.rm=T)
	n=sum(unlist(Distribution[[i]]), na.rm=T)-m
	k_list=sapply(Distribution[[i]], sum)

	Ctp_Chr_Count =sum(mapply(function(q, k) phyper(q-1, m, n, k, lower.tail=F), q_list, k_list)[H1_Events_Enrichment[[i]]] <0.4)

	if(Ctp_Chr_Count){
		Chr_CTP=rownames(as.matrix(which(mapply(function(q, k) phyper(q-1, m, n, k, lower.tail=F), q_list, k_list)[H1_Events_Enrichment[[i]]]<0.4)))
		H2_Results[[Id_Flag]]= c(Filename=as.character(CNV_Type2_List[[i]][1,2]), Chr_id= Chr_CTP)
		Id_Flag= Id_Flag + 1 
	}
} 
H2_Results = lapply(H2_Results, function(x) x[x!='X']) 
Samples_wd_Chtp=which(!is.na(match(names(CNV_Type2_List), sapply(H2_Results, function(x) x[1]))))

############################################################################################################################################################################################
#Distance Calculation
############################################################################################################################################################################################
source('/cbcb/project2-scratch/sanju/Chromotrypsis/3.Tools/BPC_Test.R')
Location_Dist = mapply(function(x, y) Format_BPC_Events(CNV_Type2_List[[x]], as.numeric(y[-1])), Samples_wd_Chtp, H2_Results)

############################################################################################################################################################################################
#Intermediate Observation:: A chtp must have More than 2 events
############################################################################################################################################################################################
##Updated based on events counts:: A chtp event must have more than one CNV change
Samples_wd_multiple_CHTP_Events = lapply(Location_Dist, function(y) sapply(y, function(x) nrow(x[2][[1]])>1) )
Samples_wd_Chtp_updated=Samples_wd_Chtp[sapply(Samples_wd_multiple_CHTP_Events, sum)>0]

Location_Dist_updated = mapply(function(x, y) x[y], Location_Dist, Samples_wd_multiple_CHTP_Events)
Location_Dist_updated= Location_Dist_updated[sapply(Location_Dist_updated, length) > 0]

Results_updated = mapply( function(x, y) c(x[1], x[-1][y]) , H2_Results, Samples_wd_multiple_CHTP_Events )
Results_updated = Results_updated[sapply(Results_updated, function(x) length(x)>1 )]

############################################################################################################################################################################################
#Hypothesis 3:: Breakpoint CLustering:: Breakpoints should be closer than background distance.
############################################################################################################################################################################################

H3_Results = lapply(Location_Dist_updated, function(x) Breakpoint_Cluster(x))   
H3_True_Samples = sapply(H3_Results, function(x) x[1,]<0.2 & sapply(x[2,], length)>2)

Results_updated_H3 = mapply(function(x, y) c(x[1], x[-1][y]) , Results_updated,  H3_True_Samples)
Results_updated_H3 =Results_updated_H3[sapply(Results_updated_H3, length)>1] 
length(Results_updated_H3)

#CHromosome enrichmetn
sort(table(unlist(sapply(Results_updated_H3, function(x) x[-1]))))
saveRDS(Results_updated_H3, '/cbcb/project2-scratch/sanju/Chromotrypsis/4.Results/Samples_1.RDS')
#It's the possibility of having a dream come true that makes life interesting.	
