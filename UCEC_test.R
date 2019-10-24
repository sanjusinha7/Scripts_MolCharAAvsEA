##TEst UCEC
cancer_type='UCEC'
ucec=readRDS('/Users/sinhas8/Project_Chromotrypsis/UCEC_infAncestry.RDS')
tcga=read.csv('/Users/sinhas8/Project_Chromotrypsis/Results_New/Nov_28/Supp_Table3.csv')
ucec_tcga=tcga[match(ucec$sample, tcga$info_tcga.Sample_Name),]

wilcox.test(ucec_tcga$CNV_Burden ~ ucec$inferred_ancestry, alternative='g')
summary(lm(ucec_tcga$CNV_Burden ~ ucec$inferred_ancestry+ucec$age, alternative='g'))
