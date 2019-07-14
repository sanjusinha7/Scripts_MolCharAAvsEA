################################################################
##Preapring LUSC race specific GISTIC analysis
################################################################
setwd('/Users/sinhas8/Project_Chromotrypsis/')
LUSC=read.csv('2.Data/Seg_TCGA/LUSC/LUSC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt',
              sep='\t')
LUAD=read.csv('2.Data/Seg_TCGA/LUAD/LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt', 
              sep='\t')
tcga=read.csv('Results_New/Nov_28/Supp_Table3.csv')
################################################################
##Pre-processing
################################################################
LUSC=LUSC[as.numeric(substring(sapply(LUSC$Sample, function(x) strsplit(as.character(x), '-')[[1]][4]), 1, 2))==1,]
LUAD=LUAD[as.numeric(substring(sapply(LUAD$Sample, function(x) strsplit(as.character(x), '-')[[1]][4]), 1, 2))==1,]

# table(as.numeric(substring(sapply(LUAD$Sample, function(x) strsplit(as.character(x), '-')[[1]][5]), 1, 2)))
# LUAD_list=split(LUAD, LUAD$Sample)
# dim(LUAD)
#  dim(do.call(rbind, lapply(LUAD_list, function(x) 
#   x[substring(sapply(x$Sample, function(x) strsplit(as.character(x), '-')[[1]][5]), 1, 2)==
#     levels(factor(substring(sapply(x$Sample, function(x) strsplit(as.character(x), '-')[[1]][5]), 1, 2)))[1], ] 
#   )))
# 
# table(unlist(sapply(LUAD_list, function(x) 
#   length(levels(factor(substring(sapply(x$Sample, 
#                                         function(x) 
#                                           strsplit(as.character(x), '-')[[1]][5]), 1, 2))))) 
# ))
# 
# table(as.numeric(substring(sapply(LUAD$Sample, function(x) strsplit(as.character(x), '-')[[1]][5]), 1, 2)))




LUSC$Sample=apply(sapply(LUSC$Sample, function(x) strsplit(as.character(x), '-')[[1]][1:3]), 
      2, 
      function(x) paste0(x, collapse = '-'))
LUAD$Sample=apply(sapply(LUAD$Sample, function(x) strsplit(as.character(x), '-')[[1]][1:3]), 
                  2, 
                  function(x) paste0(x, collapse = '-'))
LUSC_EA=LUSC[!is.na(match(LUSC$Sample,
      tcga$info_tcga.Sample_Name[tcga$info_tcga.hist=='LUSC' & tcga$info_tcga.race=='WHITE'])),]
LUSC_AA=LUSC[!is.na(match(LUSC$Sample,
                    tcga$info_tcga.Sample_Name[tcga$info_tcga.hist=='LUSC' & tcga$info_tcga.race=='AA'])),]
LUAD_EA=LUAD[!is.na(match(LUAD$Sample,
                          tcga$info_tcga.Sample_Name[tcga$info_tcga.hist=='LUAD' & tcga$info_tcga.race=='WHITE'])),]
LUAD_AA=LUAD[!is.na(match(LUAD$Sample,
                          tcga$info_tcga.Sample_Name[tcga$info_tcga.hist=='LUAD' & tcga$info_tcga.race=='AA'])),]
################################################################
##Save race-specific files
################################################################
dim(LUSC_AA)
write.table(LUSC_EA,file='/Users/sinhas8/LUSC_EA.txt', quote = F, row.names = F, sep='\t')
write.table(LUSC_AA,file='/Users/sinhas8/LUSC_AA.txt', quote = F, row.names = F, sep='\t')
write.table(LUAD_EA,file='/Users/sinhas8/LUAD_EA.txt', quote = F, row.names = F, sep='\t')
write.table(LUAD_AA,file='/Users/sinhas8/LUAD_AA.txt', quote = F, row.names = F, sep='\t')

cd ..
cd sinhas8
lcd ..
put LUSC_EA.txt
put LUSC_AA.txt
put LUAD_EA.txt
put LUAD_AA.txt

