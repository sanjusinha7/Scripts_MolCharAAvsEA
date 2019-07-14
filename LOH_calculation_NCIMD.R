##This is nor completed neither correct. 
temp=read.csv('/home/sinhas8/Downloads/From_OncoScan_COnsole/Result_Try1.segment.txt', sep='\t')
temp_LOH=temp[temp$Type == 'LOH',]

temp_filtered=temp_LOH[,c(12, 2, 3, 4, 5, 7)]
temp_filtered = temp_filtered[!(temp_filtered$Chromosome == 24 | temp_filtered$Chromosome ==25),]
temp_filtered=temp_filtered[temp_filtered$StartPosition<temp_filtered$StopPosition,]
temp_filtered=temp_filtered[is.finite(temp_filtered$State),]

#write.table(temp_filtered, '/home/sinhas8/Downloads/Segments_wdout_LOH_sex_try1_17th_complete.txt', quote=F, row.names=F, sep='\t', col.names=F)

demo=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/demo_and_clones.csv')
temp_filtered_AA=temp_filtered[unlist(sapply(demo$acc[demo$race=='AA'], function(x) grep(x, temp_filtered$FileName))),]
temp_filtered_EA=temp_filtered[unlist(sapply(demo$acc[demo$race=='EA'], function(x) grep(x, temp_filtered$FileName))),]

#write.table(temp_filtered_AA, '/home/sinhas8/Downloads/Segments_wdout_LOH_sex_try1_17th_complete_AA.txt', quote=F, row.names=F, sep='\t', col.names=F)
#write.table(temp_filtered_EA, '/home/sinhas8/Downloads/Segments_wdout_LOH_sex_try1_17th_complete_EA.txt', quote=F, row.names=F, sep='\t', col.names=F)



###LOH in NCIMD dataset
both=read.table('/home/sinhas8/Downloads/From_OncoScan_COnsole/Both_Phases_combined.txt')
colnames(both)=colnames(temp_filtered)
#both$State=log(both$State/2,2)
both$State[both$State == -Inf] = -5

demo_AA       = demo[demo$race=='AA',]
both_AA       = both[unlist(sapply(demo_AA$acc, function(x) grep(x, both$FileName))),]
both_adeno_AA = both[unlist(sapply(demo_AA$acc[demo_AA$composite_histo=='adeno'], function(x) grep(x, 
                both$FileName))),]
both_sq_AA    = both[unlist(sapply(demo_AA$acc[demo_AA$composite_histo=='sq'], function(x) grep(x, 
                both$FileName))),]

demo_EA      = demo[demo$race=='EA',]
both_EA      = both[unlist(sapply(demo_EA$acc, function(x) grep(x, both$FileName))),]
both_adeno_EA= both[unlist(sapply(demo_EA$acc[demo_EA$composite_histo=='adeno'], function(x) grep(x, 
               both$FileName))),]
both_sq_EA=both[unlist(sapply(demo_EA$acc[demo_EA$composite_histo=='sq'], function(x) grep(x, both$FileName))),]

#write.table(both_AA, '/home/sinhas8/Downloads/Both_AA.txt', quote=F, row.names=F, sep='\t', col.names=F)
#write.table(both_EA, '/home/sinhas8/Downloads/Both_EA.txt', quote=F, row.names=F, sep='\t', col.names=F)

#write.table(both_adeno_AA, '/home/sinhas8/Downloads/Both_adeno_AA.txt', quote=F, row.names=F, sep='\t', col.names=F)
#write.table(both_adeno_EA, '/home/sinhas8/Downloads/Both_adeno_EA.txt', quote=F, row.names=F, sep='\t', col.names=F)
#write.table(both_sq_AA, '/home/sinhas8/Downloads/Both_sq_AA.txt', quote=F, row.names=F, sep='\t', col.names=F)
#write.table(both_sq_EA, '/home/sinhas8/Downloads/Both_sq_EA.txt', quote=F, row.names=F, sep='\t', col.names=F)



##Calculating CNV burden 
#both_adeno_AA=both_adeno_AA[both_adeno_AA$State,]

CNV_chromosome_burden<-function(both){
	both_byfile=split(both,both$FileName)
	both_byfileandchr=sapply(both_byfile, function(x) split(x, x$Chromosome))
	new_CNV_burden=sapply(both_byfileandchr, function(x) sum(sapply(x, function(y) 
                       sum(y$StopPosition[y$State!= 2]-y$StartPosition[y$State!= 2])))/(3.3*10^9) )
	new_CNV_burden
}

GI=CNV_chromosome_burden(both)
