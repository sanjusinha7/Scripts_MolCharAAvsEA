#The aim is to write a frequency of alteration of a significant recurrent region.
#########################
##Step 0: Define Functions
#########################
###Functions required for the above job
writing_freq_file<-function(folder_location){
  setwd(folder_location)
  Filename_all='all_lesions.conf_90.txt'
  Filename_Del='del_genes.conf_90.txt'
  Filename_Amp='amp_genes.conf_90.txt'
  
  mat=read.csv(Filename_all, sep='\t')
  mat_CNV=mat[mat$Amplitude.Threshold==names(table(mat$Amplitude.Threshold)[3]),]
  mat_CNV=mat_CNV[, -(7:9)]
  sample_count=apply(mat_CNV[,-(2:6)], 1, function(x) ifelse(as.logical(grepl('Amp',x[1])), sum(as.numeric(x[-1])>0.1, na.rm=T), sum(as.numeric(x[-1])< -0.1, na.rm=T)) )
  mat_CNV$Frequency=(sample_count/(ncol(mat_CNV[,-(1:6)])-1))	
  new_mat=data.frame(mat_CNV[,1:6], Frequency=mat_CNV$Frequency)
  
  ##Adding Genes from another file
  new_mat_Del=new_mat[grep('Del',new_mat$Unique.Name),]
  new_mat_Amp=new_mat[grep('Amp',new_mat$Unique.Name),]
  
  Genes_Del=make_GeneList(Filename_Del)	
  Genes_Amp=make_GeneList(Filename_Amp)	
  
  new_mat_Del$Genes=sapply(gsub(' ', '',as.character(new_mat_Del$Descriptor)), function(x)  paste(Genes_Del$GeneName[grep(x, Genes_Del$cytoband)], collapse='| ') )
  new_mat_Amp$Genes=sapply(gsub(' ', '',as.character(new_mat_Amp$Descriptor)), function(x)  paste(Genes_Amp$GeneName[grep(x, Genes_Amp$cytoband)], collapse='| ') )
  
  rbind(new_mat_Del, new_mat_Amp)
}
##Genelist from a Amplification or Deletion output of GISTIC
make_GeneList<-function(GISTIC_output_File_id){
  Genes=read.csv(GISTIC_output_File_id, fill=T, sep='\t', header=F)
  cytoband=Genes[1,-1]
  tt=Genes[-c(1:4),-1]
  GeneScore = mapply(function(a, b) data.frame(GeneName=unique(a[a!='']),cytoband=b) ,tt, cytoband, SIMPLIFY=F)
  do.call(rbind, GeneScore)
}

#########################
##Step 0: Calling functions
#########################
##For NCIMD
adeno_EA=writing_freq_file('/home/sinhas8/Both_adeno_EA/')
sq_EA=writing_freq_file('/home/sinhas8/Both_sq_EA/')
adeno_AA=writing_freq_file('/home/sinhas8/Both_adeno_AA/')
sq_AA=writing_freq_file('/home/sinhas8/Both_sq_AA/')
##For TCGA
LUAD_AA=writing_freq_file('/Users/sinhas8/Project_Chromotrypsis/4.Results/gistic_TCGA/LUAD_AA/')
LUAD_EA=writing_freq_file('/Users/sinhas8/Project_Chromotrypsis/4.Results/gistic_TCGA/LUAD_EA/')
LUSC_AA=writing_freq_file('/Users/sinhas8/Project_Chromotrypsis/4.Results/gistic_TCGA/LUSC_AA/')
LUSC_EA=writing_freq_file('/Users/sinhas8/Project_Chromotrypsis/4.Results/gistic_TCGA/LUSC_EA/')

setwd('/Users/sinhas8/Project_Chromotrypsis/4.Results/gistic_TCGA/')
list.files()
write.csv(LUAD_EA, 'LUAD_EA/freq_of_cytobands_wd_Genes.csv')
write.csv(LUAD_AA, 'LUAD_AA/freq_of_cytobands_wd_Genes.csv')
write.csv(LUSC_EA, 'LUSC_EA/freq_of_cytobands_wd_Genes.csv')
write.csv(LUSC_AA, 'LUSC_AA/freq_of_cytobands_wd_Genes.csv')

