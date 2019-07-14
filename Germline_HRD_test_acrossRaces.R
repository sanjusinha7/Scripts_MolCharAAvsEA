###########################################################################################################################################################################
#Get Files needed
###########################################################################################################################################################################
DDR=read.csv('/Users/sinhas8/Downloads/DDR_Genes.csv')
colnames(DDR$Gene.Symbol); colnames(DDR)
DDR_Geneset_List=apply(DDR[,c(11:29)], 2, function(x) unlist(x[x!='']) )
HR_Genes=unlist(DDR[,15][DDR[,15]!=''])
core_HRGenes=unlist(DDR[,25][DDR[,25]!=''])
gtex=read.csv('/Users/sinhas8/Downloads/Gtex_Mutation.csv')
gtex=gtex[gtex$Race==2 |gtex$Race==3,] 
gtex$Race=factor(as.character(gtex$Race), labels = c('AA','EA'))
gtex$Matched.Normal.Sample.Name=factor(as.character(gtex$Matched.Normal.Sample.Name))
gtex$Tissue=factor(as.character(gtex$Tissue))
levels(gtex$Tissue)
###########################################################################################################################################################################
#Number of variants across race
###########################################################################################################################################################################
#c(head(sort(table(gtex$Hugo.Symbol), decreasing = T)), table(gtex$Hugo.Symbol)['PARP1'])
#gtex[gtex$Hugo.Symbol=='PARP1',]
gtex_l1=split(gtex, gtex$Matched.Normal.Sample.Name)
length(gtex_l1)
race=sapply(gtex_l1, function(x) x$Race[1] )
length(race)
gtex_l2=lapply(gtex_l1, function(x) split(x, x$Tissue))
dim(gtex)

Test_TotalVariants<-function(geneset_OI=DDR_Geneset_List[[1]]){
  variants_tiss_Specific = sapply(1:length(gtex_l2[[1]]), function(x) sapply(gtex_l2, function(y)  length(y[[x]]$Hugo.Symbol) ) )
  colnames(variants_tiss_Specific)=levels(gtex$Tissue)
  test1=rbind(apply(variants_tiss_Specific, 2, function(x) 
    unlist(fisher.test(aggregate(x~race,data.frame(x,race), function(x) c(sum(x), sum(x==0)))[,-1])[c('p.value','estimate')]) ),
    apply(variants_tiss_Specific, 2, function(x) c(wilcox.test(x~factor(race))$p.value, aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[,2])),
    apply(variants_tiss_Specific, 2, function(x) 
      aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[1,2]-aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[2,2] ),
    apply(variants_tiss_Specific, 2, function(x) c(aggregate(x~race,data.frame(x,race), function(x)sum(x>0) )[,2]))
  )
  rownames(test1)=c('fisher_p','odds','wicox_p','AA_Prop','EA_Prop','Prop_Ratio_AAbyEA','AA_Total_Samples','EA_Total_Samples')
  test1
}
Test_Variants<-function(geneset_OI='TP53'){
  variants_tiss_Specific=sapply(1:length(gtex_l2[[1]]), function(x) sapply(gtex_l2, function(y)  length(na.omit(match(geneset_OI,y[[x]]$Hugo.Symbol)))) ) 
  colnames(variants_tiss_Specific)=levels(gtex$Tissue)
  test1=rbind(apply(variants_tiss_Specific, 2, function(x) 
    unlist(fisher.test(aggregate(x~race,data.frame(x,race), function(x) c(sum(x), sum(x==0)))[,-1])[c('p.value','estimate')]) ),
    apply(variants_tiss_Specific, 2, function(x) c(wilcox.test(x~factor(race))$p.value, aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[,2])),
    apply(variants_tiss_Specific, 2, function(x) 
      aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[1,2]-aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[2,2] ),
    apply(variants_tiss_Specific, 2, function(x) c(aggregate(x~race,data.frame(x,race), function(x) sum(x))[,2]))
  )
  rownames(test1)=c('fisher_p','odds','wicox_p','AA_Prop','EA_Prop','Prop_Ratio_AAbyEA','AA_Total_Samples','EA_Total_Samples')
  test1
}

#Total across race
TotalGermline_across_Race=Test_TotalVariants(NA)
TotalGermline_across_Race[,TotalGermline_across_Race[1,]<0.1 & TotalGermline_across_Race[6,]<0]
TotalGermline_across_Race[,TotalGermline_across_Race[1,]<0.1 & TotalGermline_across_Race[6,]>0]
#which(TotalGermline_across_Race['AA_Prop',]>TotalGermline_across_Race['EA_Prop',])

##For all DDR 
All_DDR_Genes=unique(unlist(DDR_Geneset_List))
AllDDR_across_Race=Test_Variants(geneset_OI = All_DDR_Genes)
AllDDR_across_Race[,AllDDR_across_Race[1,]<0.1 &  ]
AllDDR_across_Race[,AllDDR_across_Race[1,]<0.1 &  TotalGermline_across_Race[6,]>0]
AllDDR_across_Race[,AllDDR_across_Race[1,]<0.1 &  TotalGermline_across_Race[6,]<0]

##For DDR pathways 
DDRGermline_across_Race = sapply(1:length(DDR_Geneset_List), function(x) Test_Variants(geneset_OI = DDR_Geneset_List[[x]]), simplify = F)
names(DDRGermline_across_Race)=names(DDR_Geneset_List)
#Reduce(intersect, DDR_Geneset_List[sapply(DDRGermline_across_Race, function(x) length(which(x[1,]<0.2)) )>0])
sapply(DDRGermline_across_Race, function(x) x[,x[1,]<0.1 &x[6,]<0 ])
sapply(DDRGermline_across_Race, function(x) x[,x[1,]<0.1 &x[6,]>0 ])

##For DDR Genes
DDRGenes_across_Race= sapply(1:length(All_DDR_Genes), function(x) Test_Variants(geneset_OI = All_DDR_Genes[x]), simplify = F)
names(DDRGenes_across_Race)=All_DDR_Genes
DDRGenes_across_Race[head(names(sort(sapply(All_DDR_Genes, function(x) nrow(gtex[gtex$Hugo.Symbol==x,])) , decreasing = T)))]
sapply(DDRGenes_across_Race, function(x) x[,x[1,]<0.1 &x[6,]<0 ] )
sapply(DDRGenes_across_Race, function(x) x[,x[1,]<0.1 &x[6,]>0 ] )
#fisher.test(matrix())

#For all the genes
avana=readRDS('/Users/sinhas8/Downloads/19Q1/avana_preprocessed_common.RDS')
AllGenes_across_Race=sapply(rownames(avana), function(x) Test_Variants(x), simplify = F)
Genes_More_mutinAA=AllGenes_across_Race[sapply(names(AllGenes_across_Race), function(x) length(AllGenes_across_Race[[x]][[1]]))>0]
names(Genes_More_mutinAA)[na.omit(match(unique(unlist(DDR_Geneset_List)), names(Genes_More_mutinAA)))]

#testing for TP53
sapply(DDR_Geneset_List, function(x) which(x=='TP53') )


###Remove Esophagus and Skin
levels(gtex$Tissue)
gtex_noMelanin = gtex[!(gtex$Tissue =='Skin' | gtex$Tissue =='Esophagus' ), ]
gtex_noMelanin$Tissue=factor(as.character(gtex_noMelanin$Tissue))
gtex_noMelanin$Matched.Normal.Sample.Name = factor(as.character(gtex_noMelanin$Matched.Normal.Sample.Name))

gtex_l1_noMelanin = split(gtex_noMelanin, gtex_noMelanin$Matched.Normal.Sample.Name)
cor.test(sapply(gtex_l1_noMelanin, nrow), Age_noMelanin)
wilcox.test(sapply(gtex_l1_noMelanin, nrow) ~ race_noMelanin)
race_noMelanin = sapply(gtex_l1_noMelanin, function(x) x$Race[1] )
Age_noMelanin = sapply(gtex_l1_noMelanin, function(x) x$Age[1] )

gtex_l2_noMelanin=lapply(gtex_l1_noMelanin, function(x) split(x, x$Tissue))

# table(race_noMelanin)
# #aggregate(Hugo.Symbol~Race,gtex_noMelanin, function(x) length(x))[,2]/table(race_noMelanin)
# aggregate(Hugo.Symbol~Race,gtex_noMelanin, function(x) length(x))[,2]/table(race_noMelanin)
# aggregate(Hugo.Symbol~Race+Age, gtex_noMelanin, function(x) length(x))
race=race_noMelanin
Test_TotalVariants<-function(geneset_OI=DDR_Geneset_List[[1]]){
  variants_tiss_Specific = sapply(1:length(gtex_l2_noMelanin[[1]]), function(x) sapply(gtex_l2_noMelanin, function(y)  length(y[[x]]$Hugo.Symbol) ) )
  colnames(variants_tiss_Specific)=levels(gtex_l1_noMelanin[[1]]$Tissue)
  test1=rbind(apply(variants_tiss_Specific, 2, function(x) 
    unlist(fisher.test(aggregate(x~race,data.frame(x,race), function(x) c(sum(x), sum(x==0)))[,-1])[c('p.value','estimate')]) ),
    apply(variants_tiss_Specific, 2, function(x) c(wilcox.test(x~factor(race))$p.value, aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[,2])),
    apply(variants_tiss_Specific, 2, function(x) 
      aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[1,2]-aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[2,2] ),
    apply(variants_tiss_Specific, 2, function(x) c(aggregate(x~race,data.frame(x,race), function(x)sum(x>0) )[,2]))
  )
  rownames(test1)=c('fisher_p','odds','wicox_p','AA_Prop','EA_Prop','Prop_Ratio_AAbyEA','AA_Total_Samples','EA_Total_Samples')
  test1
}
Test_Variants<-function(geneset_OI='TP53'){
  variants_tiss_Specific=sapply(1:length(gtex_l2_noMelanin[[1]]), function(x) sapply(gtex_l2_noMelanin, function(y)  length(na.omit(match(geneset_OI,y[[x]]$Hugo.Symbol)))) ) 
  colnames(variants_tiss_Specific)=levels(gtex_l1_noMelanin[[1]]$Tissue)
  test1=rbind(apply(variants_tiss_Specific, 2, function(x) 
    unlist(fisher.test(aggregate(x~race,data.frame(x,race), function(x) c(sum(x), sum(x==0)))[,-1])[c('p.value','estimate')]) ),
    apply(variants_tiss_Specific, 2, function(x) c(wilcox.test(x~factor(race))$p.value, aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[,2])),
    apply(variants_tiss_Specific, 2, function(x) 
      aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[1,2]-aggregate(x~race,data.frame(x,race), function(x)sum(x)/length(x))[2,2] ),
    apply(variants_tiss_Specific, 2, function(x) c(aggregate(x~race,data.frame(x,race), function(x)sum(x) )[,2]))
  )
  rownames(test1)=c('fisher_p','odds','wicox_p','AA_Prop','EA_Prop','Prop_Ratio_AAbyEA','AA_Total_Samples','EA_Total_Samples')
  test1
}
TotalGermline_across_Race=Test_TotalVariants(NA)
TotalGermline_across_Race[,TotalGermline_across_Race[1,]<0.1 & TotalGermline_across_Race[6,]<0]
TotalGermline_across_Race[,TotalGermline_across_Race[1,]<0.1 & TotalGermline_across_Race[6,]>0]

All_DDR_Genes=unique(unlist(DDR_Geneset_List))
AllDDR_across_Race=Test_Variants(geneset_OI = All_DDR_Genes)
AllDDR_across_Race[,AllDDR_across_Race[1,]<0.1]
AllDDR_across_Race[,AllDDR_across_Race[1,]<0.2 &  TotalGermline_across_Race[6,]<0]
AllDDR_across_Race[,AllDDR_across_Race[1,]<0.1 &  TotalGermline_across_Race[6,]>0]
