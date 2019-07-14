#CHAS
Phase1=read.csv('/home/sinhas8/Downloads/Segment_CHAS_Phase1.segment.txt', sep='\t')
Phase2=read.csv('/home/sinhas8/Downloads/P2_Initial_Data/CHAS_Segements_Try1.csv.txt', sep='\t')
Phase1=Phase1[,c('Use.In.Export' , 'File' , 'CN.State' , 'Type' , 'Chromosome' , 'Full.Location')]
Phase2=Phase2[,c('Use.In.Export' , 'File' , 'CN.State' , 'Type' , 'Chromosome' , 'Full.Location')]

#OncoScan
Phase1=read.csv('/home/sinhas8/Downloads/From_OncoScan_COnsole/Segment_OncoScan_Phase1.segment.txt', sep='\t')
Phase2=read.csv('/home/sinhas8/Downloads/From_OncoScan_COnsole/Segment_OncoScan_Phase2.segment.txt', sep='\t')

both=rbind(Phase1, Phase2)
both=both[both$Type=='TotalCN',]
colnames(both)
mat=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/mat_Nov09.csv')

Normal_files=levels(both$File)[is.na(match(levels(both$File), mat$File))]
both_filtered=both[-which(!is.na(match(both$File, Normal_files))),]

both_list=split(both_filtered, both_filtered$File)
both_list=both_list[-which(sapply(both_list, nrow)==0)]

LOH=data.frame(sapply(both_list, function(x)  Calculate_LOH(x)))
colnames(LOH)=c('Count_LOH')
mat=

mat$LOH=LOH

Calculate_LOH<-function(cnv){
#	Count_LOH=
	sum(cnv$Type=='LOH')
#	LOH_Prop=sum(cnv$StopPosition[cnv$Type=='LOH']-cnv$StartPosition[cnv$Type=='LOH'])/(3.3*10^9)
#	c(Count_LOH, LOH_Prop)
}


###Plotting LOH across our cohort 
tiff('/home/sinhas8/Projects/Project_Chromotrypsis/NCIMD_HRD_byLOH.tiff')
theme_set(theme_gray(base_size = 18))
ggplot(mat, aes(y=LOH, x=hist, fill=race))+geom_boxplot()+
	stat_compare_means(label = "p.signif", method = "wilcox" 
,method.args = list(alternative = "l")
)+
	labs(x='Histology', y='LOH Segments Counts')
dev.off()



###############
#CNV burden
###############
both$FileName=as.character(both$FileName)

CNV_chromosome_burden<-function(both){
	both_byfile=split(both,both$FileName)
	both_byfileandchr=lapply(both_byfile, function(x) split(x, x$Chromosome))
	CNV_burden=sapply(both_byfileandchr, function(x) sum(sapply(x, function(y) 
                       sum(y$StopPosition[y$State!= 2]-y$StartPosition[y$State!= 2])))/(3.3*10^9) )
	Loss_burden=sapply(both_byfileandchr, function(x) sum(sapply(x, function(y) 
                       sum(y$StopPosition[y$State< 2]-y$StartPosition[y$State< 2])))/(3.3*10^9) )
	Gain_burden=sapply(both_byfileandchr, function(x) sum(sapply(x, function(y) 
                       sum(y$StopPosition[y$State> 2]-y$StartPosition[y$State> 2])))/(3.3*10^9) )

	data.frame(CNV_burden, Gain_burden, Loss_burden)
}

GI=CNV_chromosome_burden(both)


