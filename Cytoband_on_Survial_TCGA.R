#Cytoband Matrix
#########################
##Step 0: Libraries
#########################
require(survival); require(survminer)
'%!in%' <- function(x,y)!('%in%'(x,y))
#########################
##Step 0: Load Files
#########################
setwd('/Users/sinhas8/Project_Chromotrypsis/')
mat=read.csv('2.Data/Corrected_NCIMD_HRD_by_LOH_and_GI.csv')
mat_updated=read.csv('2.Data/demo_updated.csv')
cytobands=read.csv('2.Data/cytobandv1', sep='\t')
#########################
##Step 0: pre-processing
#########################
mat_updated=mat_updated[na.omit(match(mat$acc, mat_updated$acc)),]
mat=mat[!is.na(match(mat$acc, mat_updated$acc)),]
mat_updated=mat_updated[,c(1, 4, 6)]
mat_updated$hist=mat$hist; mat_updated$race=mat$race
gsub('chr','',cytobands$X.chrom)
cytobands$Chr=as.numeric(gsub('chr','',sapply(cytobands$X.chrom, function(x) 
  strsplit(as.character(x), '_')[[1]][1])))
cytobands=na.omit(cytobands)
cytobands$loc=sapply(seq(nrow(cytobands)), function(x)
  paste(cytobands$Chr[x], cytobands$name[x], sep=''))
#########################
##Step 0: Define Functions
#########################
return_peaks<-function(hist='LUSC', race='AA'){
  qq=read.csv(paste('4.Results/gistic_TCGA/',
                    hist, '_',
                    race, '/',
                    'all_lesions.conf_90.txt',
                    sep=''), sep='\t')
  
  qq_FF=qq[qq$Amplitude.Threshold!=levels(qq$Amplitude.Threshold)[3],]
  qq_FF=qq_FF[,-ncol(qq_FF)]
  proc_CNV=qq_FF
  rownames(proc_CNV)=paste0(proc_CNV$Unique.Name, '_',
                            proc_CNV$Descriptor)
  proc_CNV=proc_CNV[,-c(1:9)]
  proc_CNV=t(proc_CNV)

  # surv = mat[na.omit(match(substring(rownames(proc_CNV), 2,), 
  #                          gsub('-','.',mat$File) )) ,
  #            c('lungcancer_death_all_years','survival', 'stage_trimmed')]
  proc_CNV=proc_CNV[,apply(proc_CNV, 2, function(x) length(table(x)) )>1]
  
  df=proc_CNV
  head(df)

  Region=gsub(' ','',sapply(colnames(df), 
                            function(x) strsplit(x,'_')[[1]][2]))
  df_peaks=data.frame( Frequency=apply(proc_CNV, 2, function(x) sum(x>0)/length(x)),
                      Type=sapply(colnames(df), 
                                  function(x) ifelse(grepl('Del',x), 'Del', 'Amp')),
                      Region)
  #calculate their Frequency
  df_peaks
}
Filter2_AA_specific_recurrent<-function(ForAA, ForEA){
  AA_Peaks=rownames(ForAA)
  EA_Peaks=rownames(ForEA)
  AA_Peaks=gsub(' ','',sapply(AA_Peaks, function(x) strsplit(x, '_')[[1]][2]))
  EA_Peaks=gsub(' ','',sapply(EA_Peaks, function(x) strsplit(x, '_')[[1]][2]))
  AA_Peaks_wdLocation=cbind(AA_Peaks, cytobands[match(AA_Peaks, cytobands$loc),])
  EA_Peaks_wdLocation=cbind(EA_Peaks, cytobands[match(EA_Peaks, cytobands$loc),])
  K=1
  ForAA[rownames(AA_Peaks_wdLocation[which(!as.logical(sapply(1:nrow(AA_Peaks_wdLocation), function(K)
    tenpercentRegion_test(EA_Peaks_wdLocation, AA_Peaks_wdLocation[K,])))),]),]
}
tenpercentRegion_test<-function(listA=EA_Peaks_wdLocation, region=AA_Peaks_wdLocation[K,],
                                extra=0.1){
  listA$chromStart= as.numeric(listA$chromStart)
  listA$chromEnd= as.numeric(listA$chromEnd)
  test_location= data.frame(Chr=region$Chr, 
                            start=max(0, region$chromStart-(extra*region$chromStart)),
                            end=region$chromEnd+(extra*region$chromEnd),
                            Type=rownames(region))
  listA=listA[listA$Chr==test_location$Chr & 
                (grepl('Del',rownames(listA))==grepl('Del',test_location$Type)) ,]
  sum(listA[,3]>test_location$start & listA[,3]<test_location$end | 
        listA[,4]>test_location$start & listA[,4]<test_location$end) 
}
Peaks_Freq_TCGA<-function(hist='LUSC'){
  ##Calculate their Frequency
  Freq=read.csv(paste('4.Results/gistic_TCGA/Freq_',
                      hist, sep=''))
  Freq$cyto_AA=gsub(' ','',sapply(Freq$cyto_AA, function(x) strsplit(as.character(x),'\\|')[[1]][1]))
  Freq$cyto_EA=gsub(' ','',sapply(Freq$cyto_EA, function(x) strsplit(as.character(x),'\\|')[[1]][1]))

  CNV_AA=read.csv(paste('4.Results/gistic_TCGA/',Cancer_Type, '_AA','/', 'all_data_by_genes.txt', sep=''), sep='\t', row.names=1)[,-(1:2)]
  CNV_AA=CNV_AA[-which(duplicated(sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]))),]	
  rownames(CNV_AA)=sapply(rownames(CNV_AA), function(x) strsplit(x, '\\|')[[1]][1]) 
  
  CNV_hist=cbind(CNV_EA, CNV_AA)
  #CNV_hist=CNV[,na.omit(match(mat_updated$acc[mat_updated$hist==hist], colnames(CNV)))]
  Freq_AAPeaks=t(sapply(split(Freq[,c(6,7,8,9,10,11)], Freq$cyto_AA),
                        function(x) apply(x[,-(1:2)], 2, median) ))
  Freq_EAPeaks=t(sapply(split(Freq[,c(6,7,8,9,10,11)], Freq$cyto_EA),
                        function(x) apply(x[,-(1:2)], 2, median) ))
  #For AA
  ForAA=return_peaks(hist, race='AA')
  Freq_AAPeaks=Freq_AAPeaks[na.omit(match(as.character(ForAA$Region), rownames(Freq_AAPeaks))),]
  ForAA=ForAA[!is.na(match(as.character(ForAA$Region), rownames(Freq_AAPeaks))),]
  ForAA=cbind(ForAA, Freq_AAPeaks)
  ForEA=return_peaks(hist, race='EA')
  Freq_EAPeaks=Freq_EAPeaks[na.omit(match(as.character(ForEA$Region), rownames(Freq_EAPeaks))),]
  ForEA=ForEA[!is.na(match(as.character(ForEA$Region), rownames(Freq_EAPeaks))),]
  ForEA=cbind(ForEA, Freq_EAPeaks)
  df_AAPeaks=sapply(split(Freq$X, Freq$cyto_AA),
                    function(x) apply(CNV_hist[as.character(x),], 2, median) )
  df_AA=NA
  df_EAPeaks=sapply(split(Freq$X, Freq$cyto_EA),
                    function(x) apply(CNV_hist[as.character(x),], 2, median) )
  df_EA=NA
  list(ForAA, ForEA, df_AA, df_EA)  
}


All_Steps<-function(hist='LUAD', K=1, whether_leg='none', return_plot=FALSE, numPlots=2){
  ForAA=return_peaks(hist, race='AA')
  ForEA=return_peaks(hist, race='EA')
  #Filter 1: at least 5% of the population
  ForAA=ForAA[ForAA$Frequency>0.1,]

  #Filter 2: AA-specific recurrent
  ForAA_specific=Filter2_AA_specific_recurrent(ForAA, ForEA)
  #Preprocessing
  ForAA_wdFreq=Peaks_Freq_TCGA(hist)[[1]][rownames(ForAA_specific),]
  ForAA_Amp=ForAA_wdFreq[grep('Amp',rownames(ForAA_wdFreq)),c(1:3, 4,5)]
  colnames(ForAA_Amp)[4:5]=c('Freq_in_EA', 'Freq_in_AA')
  ForAA_Del=ForAA_wdFreq[grep('Del',rownames(ForAA_wdFreq)),c(1:3, 6,7)]
  colnames(ForAA_Del)[4:5]=c('Freq_in_EA', 'Freq_in_AA')
  ForAA=rbind(ForAA_Amp, ForAA_Del)

  #Filter 3: K times more in AA
  ForAA=ForAA[ForAA$Freq_in_AA>(ForAA$Freq_in_EA*2),]
  ForAA=ForAA[order(ForAA$Freq_in_AA-ForAA$Freq_in_EA, decreasing = T),]
  ForAA=ForAA[order(ForAA$Freq_in_AA/ForAA$Freq_in_EA, decreasing = T),]
    return(ForAA)
  # }
}
#########################
##Step 0: Call Functions
#########################
Panel2=All_Steps(hist='LUSC', 
                 K=1,
                 'top', 
                 return_plot=TRUE)
leg=get_legend(All_Steps(hist='LUSC', 1, 'top', return_plot=TRUE))
######################################
####Plot Freq plot
######################################
hist='LUAD'
df_OI=data.frame(All_Steps(hist, 1, 'top', return_plot=FALSE, numPlots = 10), hist)
df2plot=data.frame(Peak=rep(paste(df_OI$Type,
                                  df_OI$Region), each=2),
                   Freq=c(apply(df_OI, 1, function(x) as.numeric(x[c(4,5)]) )),
                   race=rep(c('EA', 'AA'), nrow(df_OI)),
                   Hist=rep(df_OI$hist, each=2))
colorType_Set='Set1'
Panel1=ggplotGrob(
  ggplot(df2plot, aes(y=Freq, x=race, fill=race))+
    geom_bar(stat = 'identity')+
    labs(x="Race", y = "Frequency")+
    theme_bw(base_size = 20)+
    guides(fill=FALSE)+
    scale_fill_brewer(palette=colorType_Set)+
    facet_grid(Peak~Hist, scales = 'free', switch = 'y')+
    #geom_text(data=Pvalues, aes(x=xstar, y=ystar, label=p_values), size=8)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.y = element_text(size = 30),
          strip.text.x = element_text(size = 30)
          #strip.background.y = element_blank(),
          #strip.text.y = element_blank())
    ))

tiff(paste('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Figure3Bv2_TCGA_',hist,'Comp_Supp',Sys.Date(),'.tif'), 
     width=300, height=(nrow(df2plot)/2)*250 )
plot(Panel1)
dev.off()



######################################
####Putting em together
######################################
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Figure3Bv1.tif', 
     width=700, height=700)
lay=rbind(c(NA,2),
          c(1,2), c(1,2), c(1,2), c(1,2),
          c(1,2), c(1,2), c(1,2), c(1,2), c(1,2))
grid.arrange(Panel1, 
             plot_grid(leg, Panel2, ncol=1, rel_heights = c(1/10, 9/10)), 
             layout_matrix=lay)
dev.off()
