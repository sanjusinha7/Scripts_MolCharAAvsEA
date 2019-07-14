#################
###There is no difference btw PD-1 expression btw the two population
#################
wilcox.test(
  unlist(Exp_matched['CD274',!is.na(match(colnames(Exp_matched), mat$acc[mat$race=='AA' & mat$hist=='LUSC']))]),
  unlist(Exp_matched['CD274',!is.na(match(colnames(Exp_matched), mat$acc[mat$race=='EA' & mat$hist=='LUSC']))]),
  alternative = 'g'
)

