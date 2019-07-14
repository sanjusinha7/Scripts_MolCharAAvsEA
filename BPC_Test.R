#This script is comprised of functions to define a feature of Chromotripsis called:: Breakpoint Cluster.


#This script is comprised of functions to define a feature of Chromotripsis called:: Breakpoint Cluster in Illumina assay called:: GSA.
CNV_Types_Limit=1
BPC_Test<- function(sample){
#	CNV_Test=CNV_Test[!(CNV_Test$Chromosome=='X' | CNV_Test$Chromosome=='XY' | CNV_Test$Chromosome=='Y'),]
#	Events = split(CNV_Test, (CNV_Test$CN.State==1 | CNV_Test$CN.State==3))
	Chr2Test=names(sample)
	Distance= sapply(sample, function(y) Event_Distance(y))

	CHTP_characteristics = sapply(Chr2Test, function(x) Iterative_Wilcox_Test(
                               list(Dist=Distance[1,][x][[1]], CNV_Type = Distance[2,][x][[1]]),
                               unlist(Distance[1,][-match(x, Chr2Test)])))

#	if(length(Events)==2){
#		Chtp_Events_Dist = lapply(split(data.frame(Start=as.numeric(Events[[2]]$Start), End=as.numeric(Events[[2]]$End)), Events[[2]]$Chromosome), function(x) Event_Distance(x) )
#		lapply(Chtp_Events_Dist, function(x) Iterative_Wilcox_Test(x, NChtp_Events_Dist))
#	}
#	else{NULL}
	print('Done')
	CHTP_characteristics
}

Event_Distance<- function(DF){
	Distance_btw_Events=sapply(1:(nrow(DF)-1), function(i) DF$start[i+1]-DF$end[i])
	Copy_Number_Types=sapply(1:nrow(DF), function(i) DF$CNV[i])
	list(Distance_btw_Events, Copy_Number_Types)
}

##

##
Iterative_Wilcox_Test<- function(a, b){
	Test= tryCatch(wilcox.test(a$Dist, b, alternative='l')$p.value, error=function(cond) { 1 })

	##Key iteration
	while(length(a$Dist)>5 & Test!=1 ) {
		Test = tryCatch(wilcox.test(a$Dist, b, alternative='l')$p.value, error=function(cond) { 1 })
		if(Test>0.1 | Test<0.1 & length(table(a$CNV_Type)) > CNV_Types_Limit){
			a = Remove_Distant_Event(a)
		}
	 	else if(Test<=0.1 & length(table(a$CNV_Type)) <= CNV_Types_Limit){
			break
		}
	}
	data.frame(P.Value=Test, Total_Events_Count= length(b), CHTP_Events_Count=length(a$Dist),
                   Median_Distance=median(a$Dist))
}

Remove_Distant_Event<- function(a){
	First_Event= 1
	Last_Event= length(a$Dist)
	if(a$Dist[First_Event] > a$Dist[Last_Event]){
		a$Dist=a$Dist[-First_Event]
		a$CNV_Type=a$CNV_Type[-First_Event]
	}else{
		a$Dist=a$Dist[-Last_Event]
		a$CNV_Type=a$CNV_Type[-(Last_Event-1)]
	}	
	a
}


#Control
#mapply(function(a, b) tryCatch(wilcox.test(a, b, paired= F, correct=F, exact=F, alternative='l')$p.value, error=function(cond) { 1 }) ,  Distance_CHTP_Event, Distance_N.CHTP_Event)	
#New
#mapply(function(a, b) Iterative_Wilcox_Test(a, b) ,  Distance_CHTP_Event, Distance_N.CHTP_Event)	

