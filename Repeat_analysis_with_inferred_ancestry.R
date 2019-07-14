##Repeat the analysis with inferred ancestry
##for CNV burden
wilcox.test(mat$CNV_burden[mat$hist=='LUAD']~mat$race[mat$hist=='LUAD'], alternative='g')
wilcox.test(mat$[mat$hist=='LUSC']~mat$race[mat$hist=='LUSC'], alternative='g')
wilcox.test(mat$CNV_burden[mat$hist=='LUAD']~mat$inferred_ancestry[mat$hist=='LUAD'], alternative='g')
wilcox.test(mat$CNV_burden[mat$hist=='LUSC']~mat$inferred_ancestry[mat$hist=='LUSC'], alternative='g')

####for HRD
wilcox.test(mat$HRD_by_LOH[mat$hist=='LUAD']~mat$race[mat$hist=='LUAD'], alternative='g')
wilcox.test(mat$HRD_by_LOH[mat$hist=='LUSC']~mat$race[mat$hist=='LUSC'], alternative='g')
wilcox.test(mat$HRD_by_LOH[mat$hist=='LUAD']~mat$inferred_ancestry[mat$hist=='LUAD'], alternative='g')
wilcox.test(mat$HRD_by_LOH[mat$hist=='LUSC']~mat$inferred_ancestry[mat$hist=='LUSC'], alternative='g')

##
fisher.test(cbind(table(mat$chtp_Quan[mat$hist=='LUSC' & mat$race=='EA']),
                  table(mat$chtp_Quan[mat$hist=='LUSC' & mat$race=='AA'])),
            alternative = 'g')

write.csv(mat,'/Users/sinhas8/Project_Chromotrypsis/2.Data/mat_NCIMD_latest_27thMay.csv')
