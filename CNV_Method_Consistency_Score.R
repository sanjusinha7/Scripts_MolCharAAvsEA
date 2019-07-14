#This script is to compare the two methods: GSA & Oncoscan; Where the CNV data is compared for similarity
library(bedr)
##Files needed
CNV_GSA=readRDS('/cbcb/project2-scratch/sanju/Chromotrypsis/2.Data/CG0314_GSA_Initial_standard_cnv_report.RDS')
CNV_onco=readRDS('/cbcb/project2-scratch/sanju/Chromotrypsis/2.Data/Processes_oncoscan_CN.txt')


#Removing the CNV=2 // LOH regions for seperate quality comparison
CNV_onco_excludedLOH=CNV_onco[-which(CNV_onco$Raw.CN.State==2),]
CNV_GSA_excludedLOH=CNV_GSA[-which(CNV_GSA$CN.State==2),]

#Adding whether an event is a gain or a loss.
CNV_onco_excludedLOH$Type = CNV_onco_excludedLOH$CN.State > 2
CNV_GSA_excludedLOH$Type  = CNV_GSA_excludedLOH$CN.State  > 2

##Deblinding the samples
levels(CNV_GSA_excludedLOH[,1])=tt[match(levels(CNV_GSA_excludedLOH$'Sample ID'), tt[,1]), 2]
levels(CNV_onco_excludedLOH[,1])=unlist(sapply(sapply(levels(CNV_onco_excludedLOH[,1]), function(x) strsplit(x, "\\.")[[1]][1]), function(x) strsplit(x, '_')))

#Correcting the DF to bedformat
colnames(CNV_onco_excludedLOH)[2:4]=c('chr', 'start', 'end')
colnames(CNV_GSA_excludedLOH)[2:4]=c('chr', 'start', 'end')

#Changing data types to the correct ones.
CNV_onco_excludedLOH$chr=as.character(CNV_onco_excludedLOH$chr)
CNV_GSA_excludedLOH$chr=as.character(CNV_GSA_excludedLOH$chr)

#Splitting the dataframe based on sample
CNV_GSA_ProperFormat=split(CNV_GSA_excludedLOH[,c(2:4, 8)], CNV_GSA_excludedLOH$map )
CNV_onco_ProperFormat=split(CNV_onco_excludedLOH[,c(2:4, 7)], CNV_onco_excludedLOH[,1] )

#Dividing based on CN.event type
CNV_GSA_ProperFormat = lapply(CNV_GSA_ProperFormat, function(x) split(x[,-4], x[,4] ))
CNV_onco_ProperFormat = lapply(CNV_onco_ProperFormat, function(x) split(x[,-4], x[,4] ))

#Removing normal samples
CNV_GSA_FF=CNV_GSA_ProperFormat[-c(grep('n',names(CNV_GSA_ProperFormat)), grep('N', names(CNV_GSA_ProperFormat)))]
CNV_onco_FF=CNV_onco_ProperFormat[-c(grep('n',names(CNV_onco_ProperFormat)), grep('N', names(CNV_onco_ProperFormat)))]

##Filenames 
GSA_names=sapply(sapply(names(CNV_GSA_FF), function(x) strsplit(x, 't')[[1]][1]), function(x) strsplit(x, 'T')[[1]][1])
onco_names=sapply(sapply(names(CNV_onco_FF), function(x) strsplit(x, 't')[[1]][1]), function(x) strsplit(x, 'T')[[1]][1])

CNV_onco_SF=CNV_onco_FF[match(GSA_names, onco_names)]

R= mapply(function(a, b) mapply(function(x, y) F1(x, y) , a, b), CNV_onco_SF, CNV_GSA_FF)
sapply(R, function(x) sum(unlist(x)>30) )

#**#Useful Functions**
#FOrmatting our input
F1<- function(tt1, tt2){
	tt1[,2:3] = apply(tt1[,2:3], 2, as.numeric)
	tt2[,2:3] = apply(tt2[,2:3], 2, as.numeric)
	Chr=levels(factor(c(tt1[,1], tt2[,1])))
	sapply(Chr, function(x) RO(tt1[tt1[,1]==x,2:3], tt2[tt2[,1]==x,2:3] ) )
}

##reciprocal overlap function
RO<- function(a, b){
	OV=apply(a, 1, function(x) apply(b, 1, function(y) 2*(min(x[2], y[2])- max(x[1], y[1]))  / ( (x[2]-x[1])+(y[2]-y[1]) ) ))
	OV[OV<0]=0
	OV*100
}

 ##This script is dedicated to Caffeine and Pierre Joseph Pelletier, who chemiaclly isolated it for the first time.
##Also, I commit to provide an acknowledgement to caffiene in the paper based on above script. ;) 


