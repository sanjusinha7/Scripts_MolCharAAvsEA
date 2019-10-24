##########GEne specific enrichemnt in AA
##########################################################################################
#Custom func
##########################################################################################
Enrichmentof_Germline_Deficiency_byRace<-function(Geneset_OI=as.character(unlist(il6$Genes)),cancer_type='LUAD'){
  Geneset_OI=as.character(unlist(Geneset_OI))
  totalAA=table(demo$race[demo$acronym==cancer_type])['AA']
  totalEA=table(demo$race[demo$acronym==cancer_type])['EA']
  GermlineHRD_Total=length(na.omit(match(df$HUGO_Symbol[df$cancer==cancer_type],
                                         Geneset_OI)))
  GermlineHRD_inAA=length(na.omit(match(df$HUGO_Symbol[df$race=='AA' &
                                                         df$cancer==cancer_type],
                                        Geneset_OI)))
  GermlineHRD_inEA=length(na.omit(match(df$HUGO_Symbol[df$race=='EA' & 
                                                         df$cancer==cancer_type],
                                        Geneset_OI)))
  
  pvalue=phyper(GermlineHRD_inAA-1, totalAA, totalEA, GermlineHRD_Total, lower.tail = F)
  if(pvalue==0){
    return(data.frame(NA, NA, NA))
  } else{ return(data.frame(pvalue, GermlineHRD_inAA/totalAA, GermlineHRD_inEA/totalEA ))}
} 
PanCan_Enrichmentof_Germline_Deficiency_byRace<-function(Geneset_OI=HR_Genes){
  Geneset_OI=as.character(unlist(Geneset_OI))
  totalAA=table(demo$race)['AA']
  totalEA=table(demo$race)['EA']
  GermlineHRD_Total=length(na.omit(match(df$HUGO_Symbol,
                                         Geneset_OI)))
  GermlineHRD_inAA=length(na.omit(match(df$HUGO_Symbol[df$race=='AA'],
                                        Geneset_OI)))
  GermlineHRD_inEA=length(na.omit(match(df$HUGO_Symbol[df$race=='EA'],
                                        Geneset_OI)))
  pvalue=phyper(GermlineHRD_inAA-1, 
                totalAA, totalEA,
                GermlineHRD_Total, 
                lower.tail = F)
  if(pvalue==0){
    return(NA)
  } else{
    return(data.frame(pvalue, GermlineHRD_inAA/totalAA, GermlineHRD_inEA/totalEA ))
  }
  
}

##########################################################################################
#Call Func
##########################################################################################
# CanSpec_df2write=lapply(DDR_Geneset_List, 
#                         function(x) sapply(levels(demo$acronym), function(y) 
#                           Enrichmentof_Germline_Deficiency_byRace(Geneset_OI = x,
#                                                                   cancer_type = y)))
#######For LUSC
pathogenicvariants_HRDGenes=t(sapply(DDR_Geneset_List$Homology.dependent.recombination..HDR., 
       function(x) Enrichmentof_Germline_Deficiency_byRace(Geneset_OI = x,
                                        cancer_type = 'LUSC')))
Enrichmentof_Germline_Deficiency_byRace(Geneset_OI = c('BRCA1', 'BRCA2'),
                                        cancer_type = 'LUSC')[2]/Enrichmentof_Germline_Deficiency_byRace(Geneset_OI = c('BRCA1', 'BRCA2'),
                                        cancer_type = 'LUSC')[3]

pathogenicvariants_HRDGenes=pathogenicvariants_HRDGenes[order(unlist(pathogenicvariants_HRDGenes[,1])),]
colnames(pathogenicvariants_HRDGenes)= c('Significance of Enrichment in AA', '% AA Patients with pathogenic variant', '% EA Patients with pathogenic variant')
write.csv(pathogenicvariants_HRDGenes, 'pathogenicvariants_HRDGenes.csv')

###For PAn-Cancer

pathogenicvariants_HRDGenes_PC=t(sapply(DDR_Geneset_List$Homology.dependent.recombination..HDR., 
                                     function(x) PanCan_Enrichmentof_Germline_Deficiency_byRace(Geneset_OI = x)))
pathogenicvariants_HRDGenes_PC=pathogenicvariants_HRDGenes_PC[order(unlist(pathogenicvariants_HRDGenes_PC[,1])),]
head(pathogenicvariants_HRDGenes_PC, 15)
colnames(pathogenicvariants_HRDGenes_PC)= c('Pan Cancer-Significance of Enrichment in AA', '% AA Patients with pathogenic variant', '% EA Patients with pathogenic variant')
write.csv(pathogenicvariants_HRDGenes_PC, 'pathogenicvariants_HRDGenes_PC.csv')
