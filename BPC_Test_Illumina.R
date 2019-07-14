#This script is comprised of functions to define a feature of Chromotripsis called:: Breakpoint Cluster in Illumina assay called:: GSA.

BPC_Test<- function(CNV_Test, Events_Possible, Events_Threshold){
	CNV_Test=CNV_Test[!(CNV_Test$Chromosome=='X' | CNV_Test$Chromosome=='XY' | CNV_Test$Chromosome=='Y'),]

	pp=lapply(rownames(table(as.character(CNV_Test$Chromosome))), function(x) split(CNV_Test, CNV_Test$Chromosome==x))

	names(pp)=rownames(table(as.character(CNV_Test$Chromosome)))



	##result_Details:: Critical for Visualization
	Result_Details=lapply(pp, function(x) Iterative_Wilcox_Test( Event_Distance(OSC_Events( tryCatch(x[[2]], error=function(cond){data.frame()}), Events_Possible ) ), Event_Distance( x[[1]])  ) )



	names(Result_Details)[sapply(Result_Details, function(x) (x[3]>Events_Threshold & x[1]<0.1) & x[2]>5 )]			
}


#Returns Oscillatory CNV events
OSC_Events<- function(Events_DF, Events_Possible){
		Events_DF[!is.na(match(Events_DF$CN.State, Events_Possible)),]
}

##Given a dataframe of Start and End, following function generates distance between events.
Event_Distance<- function(DF){
	if(nrow(DF)>1){
	sapply(1:(nrow(DF)-1), function(i) as.numeric(DF$Start[i+1]) - as.numeric(DF$End[i]) )
	}	
}

#Given a Chromosome CN Breakpoint vector, this functions find the cluster which having significantly closer breakpoints.
Iterative_Wilcox_Test<- function(a, b){
	Test= tryCatch(wilcox.test(a, b, alternative='l')$p.value, error=function(cond) { 1 })
	while(Test>0.1 & Test!=1) {
		a = Remove_Distant_Event(a)
		Test= tryCatch(wilcox.test(a, b, paired= F, correct=F, exact=F, alternative='l')$p.value, error=function(cond) { 1 })
	}
	c(P.Value=Test, Total_Events_Count= length(b), CHTP_Events_Count=length(a), Median=median(a))
}

#Given a distance vector, this one removes the first or the last events:: Whichever is great.
Remove_Distant_Event<- function(a){
	First_Event= 1
	Last_Event= length(a)
	if(a[First_Event] > a[Last_Event]){
		a=a[-First_Event]
	}
	else{
		a=a[-Last_Event]
	}	
	a
}


