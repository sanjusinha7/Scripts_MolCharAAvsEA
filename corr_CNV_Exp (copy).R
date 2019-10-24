####Script to get gene level copy number from segment level:: (TCGA level 3 from TCGA level 2 file.)

##libraries
require(CNTools)
data("geneInfo")
data("sampleData")
temp = read.csv('/cbcbhomes/sanju/Segments_OncoScan.txt', sep='\t')
temp= temp[is.finite(temp[,6]),]

colnames(temp)=colnames(sampleData)

geneInfo_trunc <- geneInfo[sample(1:nrow(geneInfo), 1000), ]
cnseg=CNSeg(temp)

### Copy number change matrix by genes
rdByGene <- getRS(cnseg, by = "gene", imput = FALSE, XY = FALSE, geneMap = geneInfo, what='mean')

##Loading TCGA data for genenames::
require(data.table)
require('CNTools')

#### Presave:: Copy number change matrix by genes
CNV=rs(readRDS('/cbcbhomes/sanju/rdByGene.RDS'))
CNV_genelist=CNV[,5]

##Expression data:: Pre-saved:: using RSEM values (** could be improved by using RPKM instead of RSEM)
Exp=fread('/cbcbhomes/sanju/RawCountFile_RSEM_genes_filtered.txt')
Exp=as.matrix(Exp)
Exp_genelist= sapply(Exp[,1], function(x) strsplit(x, '\\|')[[1]][2])

Onco_Amp=readRDS('/cbcb/project2-scratch/sanju/Chromotrypsis/2.Data/Onco_Amp_genelist.RDS')
Onco_Del=readRDS('/cbcb/project2-scratch/sanju/Chromotrypsis/2.Data/Onco_Del_genelist.RDS')

Onco_Amp_Exp=Exp[na.omit(match(Onco_Amp[,1] ,Exp_genelist)),]
Onco_Del_Exp=Exp[na.omit(match(Onco_Del[,1] ,Exp_genelist)),]
Onco_Amp_CNV=CNV[na.omit(match(Onco_Amp[,1] ,CNV_genelist)),]
Onco_Del_CNV=CNV[na.omit(match(Onco_Del[,1] ,CNV_genelist)),]
 

rownames(Onco_Amp_Exp) = sapply(Onco_Amp_Exp[,1], function(x) strsplit(x, '\\|')[[1]][2])
rownames(Onco_Del_Exp) = sapply(Onco_Del_Exp[,1], function(x) strsplit(x, '\\|')[[1]][2])

##two matrices having genes amplified CNV and Expression!!
mat1 = Onco_Amp_CNV[!is.na(match(Onco_Amp_CNV[,5], rownames(Onco_Amp_Exp))),]
mat2 = Onco_Amp_Exp[na.omit(match(Onco_Amp_CNV[,5], rownames(Onco_Amp_Exp))),]

rownames(mat1)=mat1[,5]
mat1=mat1[,-(1:5)]

##
Samples_names_mat2       = sapply(sapply(colnames(mat2), function(x) strsplit(x, 'genes.')[[1]][2]), function(x) strsplit(x, '_')[[1]][1]) 
Samples_names_mat2_Tumor = Samples_names_mat2[grep('T',Samples_names_mat2)]
Samples_names_mat2_Tumor = substring(Samples_names_mat2_Tumor, 2)
mat2_Tumor               = mat2[,grep('T',Samples_names_mat2)]

##
Samples_names_mat1       = sapply(sapply(sapply(colnames(mat1), function(x) strsplit(x, '\\.')[[1]][1]), function(x) strsplit(x, 't')[[1]][1]), function(x) strsplit(x, 'T')[[1]][1])
Samples_names_mat1_Tumor = Samples_names_mat1[which(!is.na(as.double(Samples_names_mat1)))]
mat1_Tumor               = mat1[,which(!is.na(as.double(Samples_names_mat1)))]


##Matched Matrix for Amplified genes
mat2_Tumor_matched=mat2_Tumor[,which(!is.na(match(Samples_names_mat2_Tumor, Samples_names_mat1_Tumor)))]
mat1_Tumor_matched=mat1_Tumor[,na.omit(match(Samples_names_mat2_Tumor, Samples_names_mat1_Tumor))]

saveRDS(mat2_Tumor_matched,'/cbcbhomes/sanju/mat2_Tumor_matched.RDS')
saveRDS(mat1_Tumor_matched,'/cbcbhomes/sanju/mat1_Tumor_matched.RDS')

##############***********Checkpoint1:: mat1-CNV :: mat2-Expression;************##############

##
mat2_Tumor_matched=readRDS('/cbcbhomes/sanju/mat2_Tumor_matched.RDS')
mat1_Tumor_matched=readRDS('/cbcbhomes/sanju/mat1_Tumor_matched.RDS')

Genes_Amp=rownames(mat2_Tumor_matched)[which(p.adjust(sapply(1:nrow(mat2_Tumor_matched), function(x)  cor.test(as.double(mat1_Tumor_matched[x,]), as.double(mat2_Tumor_matched[x, ]), method='spearman')$p.value), method='fdr') < 0.25)]
p.value = sapply(1:nrow(mat1_Tumor_matched), function(x)  p.adjust(cor.test(as.double(mat1_Tumor_matched[x,]), as.double(mat2_Tumor_matched[x, ]), method='spearman')$p.value, method='fdr'))

require('CNTools')
df1=data.frame(genename=rownames(mat1_Tumor_matched), p.value=p.value)
CNV=rs(readRDS('/cbcbhomes/sanju/rdByGene.RDS'))

df1_wd_chr=data.frame(df1,CNV[na.omit(match(df1[,1], CNV[,5])),1:3])
df1_wd_chr_Amplified=df1_wd_chr[order(df1_wd_chr[,2]),]
GISTIC_listed_Amplified_genes_correlated_with_expression = df1_wd_chr_Amplified[which(df1_wd_chr_Amplified[,2]<0.1),]

##saving the Amplified genes corr with Expression.
write.csv(GISTIC_listed_Amplified_genes_correlated_with_expression, '/cbcb/project2-scratch/sanju/Chromotrypsis/4.Results/GISTIC_listed_Amplified_genes_correlated_with_expression.csv', quote=F, row.names=F)


###########*********************Landscape Deletion**************************##########
Onco_Del_Exp=Exp[na.omit(match(Onco_Del[,1] ,Exp_genelist)),]
Onco_Del_CNV=CNV[na.omit(match(Onco_Del[,1] ,CNV_genelist)),]

rownames(Onco_Del_Exp) = sapply(Onco_Del_Exp[,1], function(x) strsplit(x, '\\|')[[1]][2])

mat1 = Onco_Del_CNV[!is.na(match(Onco_Del_CNV[,5], rownames(Onco_Del_Exp))),]
mat2 = Onco_Del_Exp[na.omit(match(Onco_Del_CNV[,5], rownames(Onco_Del_Exp))),]

#rownames(mat1)=mat1[,5]
mat1_rownames=mat1[,5]
mat1=mat1[,-(1:5)]

##
Samples_names_mat2       = sapply(sapply(colnames(mat2), function(x) strsplit(x, 'genes.')[[1]][2]), function(x) strsplit(x, '_')[[1]][1]) 
Samples_names_mat2_Tumor = Samples_names_mat2[grep('T',Samples_names_mat2)]
Samples_names_mat2_Tumor = substring(Samples_names_mat2_Tumor, 2)
mat2_Tumor               = mat2[,grep('T',Samples_names_mat2)]

##
Samples_names_mat1       = sapply(sapply(sapply(colnames(mat1), function(x) strsplit(x, '\\.')[[1]][1]), function(x) strsplit(x, 't')[[1]][1]), function(x) strsplit(x, 'T')[[1]][1])
Samples_names_mat1_Tumor = Samples_names_mat1[which(!is.na(as.double(Samples_names_mat1)))]
mat1_Tumor               = mat1[,which(!is.na(as.double(Samples_names_mat1)))]

##Matched Matrix for Amplified genes
mat2_Tumor_matched=mat2_Tumor[,which(!is.na(match(Samples_names_mat2_Tumor, Samples_names_mat1_Tumor)))]
mat1_Tumor_matched=mat1_Tumor[,na.omit(match(Samples_names_mat2_Tumor, Samples_names_mat1_Tumor))]

Genes_Del = rownames(mat2_Tumor_matched)[which(sapply(1:nrow(mat2_Tumor_matched), function(x)  cor.test(as.double(mat1_Tumor_matched[x,]), as.double(mat2_Tumor_matched[x, ]), method='spearman')$p.value)< 0.2)]
p.value   = sapply(1:nrow(mat1_Tumor_matched), function(x)  p.adjust(cor.test(as.double(mat1_Tumor_matched[x,]), as.double(mat2_Tumor_matched[x, ]), method='spearman')$p.value, method='fdr') )

df1=data.frame(genename=rownames(mat2_Tumor_matched), p.value=p.value)
CNV=rs(readRDS('/cbcbhomes/sanju/rdByGene.RDS'))

df1_wd_chr=data.frame(df1,CNV[na.omit(match(df1[,1], CNV[,5])),1:3])
df1_wd_chr_Deleted=df1_wd_chr[order(df1_wd_chr[,2]),]

GISTIC_listed_Deleted_genes_correlated_with_expression = df1_wd_chr_Deleted[which(df1_wd_chr_Deleted[,2]<0.1),]
write.csv(GISTIC_listed_Deleted_genes_correlated_with_expression, '/cbcb/project2-scratch/sanju/Chromotrypsis/4.Results/GISTIC_listed_Deleted_genes_correlated_with_expression.csv', quote=F, row.names=F)

