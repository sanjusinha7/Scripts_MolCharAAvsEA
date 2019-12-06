###########################
####Figure 3X-Global CNV plot
########################
require(copynumber)
setwd('/Users/sinhas8/Project_Chromotrypsis/')
df=read.csv('2.Data/Input_for_GISTIC.txt',sep=' ')
Phase1=read.csv('2.Data/GISTIC_Input_forOncoScan_Raw.probeset.txt',sep='\t')
# Phase1=fread('2.Data/GISTIC_Input_forOncoScan_Raw.probeset.txt',sep='\t')
####
Phase2=cbind(read.csv('/Users/sinhas8/Downloads/Result_probeset/Result_P1.probeset.txt', sep = '\t'),
             read.csv('/Users/sinhas8/Downloads/Result_probeset/Result_P2.probeset.txt', sep = '\t')[,-(1:3)],
             read.csv('/Users/sinhas8/Downloads/Result_probeset/Result_P3.probeset.txt', sep = '\t')[,-(1:3)],
             read.csv('/Users/sinhas8/Downloads/Result_probeset/Result_P4.probeset.txt', sep = '\t')[,-(1:3)],
             read.csv('/Users/sinhas8/Downloads/Result_probeset/Result_P5.probeset.txt', sep = '\t')[,-(1:3)],
             read.csv('/Users/sinhas8/Downloads/Result_probeset/Result_P6.probeset.txt', sep = '\t')[,-(1:3)],
             read.csv('/Users/sinhas8/Downloads/Result_probeset/Result_P7.probeset.txt', sep = '\t')[,-(1:3)],
             read.csv('/Users/sinhas8/Downloads/Result_probeset/Result_P8.probeset.txt', sep = '\t')[,-(1:3)])
both=cbind(Phase1, Phase2)
saveRDS(both, '2.Data/Results_Probeset.RDS')
both=readRDS('2.Data/Results_Probeset.RDS')
mat=read.csv('2.Data/Corrected_NCIMD_HRD_by_LOH_and_GI.csv')
mat=mat[mat$hist=='adeno' | mat$hist =='sq',]
mat$hist=as.character(mat$hist)
mat$hist=factor(mat$hist)
levels(mat$hist)=c('LUAD', 'LUSC')
########################
####Renaming Column Names
########################
both_FF=both[,c(1:3,grep('WeightedLog2Ratio',colnames(both)))]
rownames(both_FF)=both_FF[,1]; both_FF=both_FF[,-1]
both_FF=both_FF[,-grep('n\\.',colnames(both_FF))]
both_FF=both_FF[,-grep('N\\.',colnames(both_FF))]
both_FF=both_FF[,-grep('Rename.',colnames(both_FF))]
colnames(both_FF)[-(1:2)] =sapply(colnames(both_FF)[-(1:2)], function(x) strsplit(x, '\\.')[[1]][3])
colnames(both_FF)[grep('_',colnames(both_FF))[1]]= substring(colnames(both_FF)[grep('_',colnames(both_FF))[1]], 
                                                             1, 3)
colnames(both_FF)[grep('_',colnames(both_FF))]=sapply(colnames(both_FF)[grep('_',colnames(both_FF))], 
                                                      function(x) strsplit(x, '_')[[1]][2])
colnames(both_FF)[-(1:2)] =sapply(colnames(both_FF)[-(1:2)], function(x) gsub("[^0-9]", "", x)  )
#Choosing Race
both_FF=both_FF[,c(1, 2, na.omit(match(mat$acc, colnames(both_FF))))]
both_FF=both_FF[-c(which(both_FF$Chromosome=='24'), which(both_FF$Chromosome=='25')),]
mat=mat[!is.na(match(mat$acc, colnames(both_FF))),]
########
race='AA'; cancer_type='LUSC'
both_FF_AA_LUSC=both_FF[,c(1,2, which(mat$race==race & mat$hist==cancer_type)+2 )]
both_FF.res_AA_LUSC <- pcf(data=both_FF_AA_LUSC, gamma=50, verbose=FALSE, fast = T)
multiseg_AA_LUSC=multipcf(data = both_FF_AA_LUSC, gamma = 100, verbose = FALSE, fast = T)
race='EA'; cancer_type='LUSC'
both_FF_EA_LUSC=both_FF[,c(1,2, which(mat$race==race & mat$hist==cancer_type)+2 )]
both_FF.res_EA_LUSC <- pcf(data=both_FF_EA_LUSC, gamma=50, verbose=FALSE, fast=T)
multiseg_EA_LUSC=multipcf(data = both_FF_EA_LUSC, gamma = 100, verbose = FALSE, fast=T)
race='AA'; cancer_type='LUAD'
both_FF_AA_LUAD=both_FF[,c(1,2, which(mat$race==race & mat$hist==cancer_type)+2 )]
both_FF.res_AA_LUAD <- pcf(data=both_FF_AA_LUAD, gamma=50, verbose=FALSE, fast=T)
multiseg_AA_LUAD=multipcf(data = both_FF_AA_LUAD, gamma = 100, verbose = FALSE, fast=T)
race='EA'; cancer_type='LUAD'
both_FF_EA_LUAD=both_FF[,c(1,2, which(mat$race==race & mat$hist==cancer_type)+2 )]
both_FF.res_EA_LUAD <- pcf(data=both_FF_EA_LUAD, gamma=50, verbose=FALSE, fast=T)
multiseg_EA_LUAD=multipcf(data = both_FF_EA_LUAD, gamma = 100, verbose = FALSE, fast=T)


Global_CNV_plot<-function(both_FF.res, multiseg, race, cancer_type){
  nseg = nrow(multiseg)
  cormat = cor(t(multiseg[,-c(1:5)]))
  chr.from <- c();pos.from <- c();chr.to <- c();pos.to <- c();cl <- c(); Pear_cor<- c()
  thresh_Pos = 0.75; thresh_Neg = -0.5
  for (i in 1:(nseg-1)) {
    for (j in (i+1):nseg) {
      #Check if segment-correlation is larger than threshold and that the two 
      #segments are located on different chromosomes
      if ((cormat[i,j] > thresh_Pos | cormat[i,j] < thresh_Neg) && multiseg$chrom[i] != multiseg$chrom[j]) {
        chr.from = c(chr.from,multiseg$chrom[i])
        chr.to = c(chr.to,multiseg$chrom[j])
        pos.from = c(pos.from,(multiseg$start.pos[i] + multiseg$end.pos[i])/2)
        pos.to = c(pos.to,(multiseg$start.pos[j] + multiseg$end.pos[j])/2)
        Pear_cor = c(Pear_cor, cormat[i,j])
        if(cormat[i,j] > thresh_Pos){
          cl <- c(cl,1)           #class 1 for those with positive correlation
        }else{
          cl <- c(cl,2)           #class 2 for those with negative correlation
        }    
      }
    }
  }
  arcs <- data.frame(chr.from,pos.from,chr.to,pos.to,cl, Pear_Cor=Pear_cor); dim(arcs)
  arcs_Pos <- arcs[arcs$Pear_Cor>0,]
  arcs_Pos =  head(arcs_Pos[order(arcs_Pos$Pear_Cor, decreasing = T),], 50)
  arcs_Neg <- arcs[arcs$Pear_Cor<0,]
  arcs_Neg =  head(arcs_Neg[order(arcs_Neg$Pear_Cor, decreasing = T),], 50)
  arcs=rbind(arcs_Neg, arcs_Pos)
  # tiff(paste('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Circ_test1_',race,'_',cancer_type,'.tif', sep=''))
  plotCircle(segments=both_FF.res,thres.gain=0.2, thres.loss = -0.2, arcs = arcs,
             freq.colors = c("blue","limegreen"))+
    title(main = paste(race, cancer_type, sep=' '), cex.main=1.8,  line = -1.5)
  # dev.off()
  # tiff(paste('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Flat_test1_',race,'_',cancer_type,'.tif', sep=''))
  # p2=plotFreq(segments=both_FF.res,thres.gain=0.2,thres.loss=-0.2)
  # dev.off()
}
tiff(paste('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Circ_test1v2.tif'),
     height = 1200, width=1200)
par(mfrow=c(2,2))
Global_CNV_plot(both_FF.res_AA_LUAD, multiseg_AA_LUAD, race='AA', cancer_type='LUAD')
Global_CNV_plot(both_FF.res_EA_LUAD, multiseg_EA_LUAD, race='EA', cancer_type='LUAD')
Global_CNV_plot(both_FF.res_AA_LUSC, multiseg_AA_LUSC, race='AA', cancer_type='LUSC')
Global_CNV_plot(both_FF.res_EA_LUSC, multiseg_EA_LUSC, race='EA', cancer_type='LUSC')
dev.off()

##############
# Vector form
##############
setEPS()
postscript(paste('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Circ_test1v2.eps'),
           height = 16, width=16)
par(mfrow=c(2,2))
Global_CNV_plot(both_FF.res_AA_LUAD, multiseg_AA_LUAD, race='AA', cancer_type='LUAD')
Global_CNV_plot(both_FF.res_EA_LUAD, multiseg_EA_LUAD, race='EA', cancer_type='LUAD')
Global_CNV_plot(both_FF.res_AA_LUSC, multiseg_AA_LUSC, race='AA', cancer_type='LUSC')
Global_CNV_plot(both_FF.res_EA_LUSC, multiseg_EA_LUSC, race='EA', cancer_type='LUSC')
dev.off()
