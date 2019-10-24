########################################################################
##BRCAness Expression signature
########################################################################
require(statar)
process_mat<-function(Exp=Exp_AA_AD){
  rownames(Exp)=make.names(Exp[,1], unique = T)
  Exp=Exp[,-1]
  myhead(Exp)
  Exp=Exp[,-grep('T',colnames(Exp))]
  colnames(Exp)=sapply(colnames(Exp), function(x) gsub("[^0-9]", "", x))
  Exp  
}
Exp_AA_AD=process_mat(read.csv('2.Data/Exp_AA_AD.csv'))
Exp_EA_AD=process_mat(read.csv('2.Data/Exp_EA_AD.csv'))
Exp_AA_SC=process_mat(read.csv('2.Data/Exp_AA_SC.csv'))
Exp_EA_SC=process_mat(read.csv('2.Data/Exp_EA_SC.csv'))

Exp=cbind(Exp_AA_AD,
      Exp_EA_AD,
      Exp_AA_SC,
      Exp_EA_SC)

demo=rbind(data.frame(id=1:ncol(Exp_AA_AD), Hist='AD', Race='AA'),
           data.frame(id=1:ncol(Exp_EA_AD), Hist='AD', Race='EA'),
           data.frame(id=1:ncol(Exp_AA_SC), Hist='SC', Race='AA'),
           data.frame(id=1:ncol(Exp_EA_SC), Hist='SC', Race='EA'))

signature=read.csv('/Users/sinhas8/Downloads/BRCAness.csv')
signature=as.character(signature$GeneName[!signature$GeneName==''])
signature_id=na.omit(match(signature, rownames(Exp)))

Exp_hrd=t(apply(Exp[signature_id,], 1, function(x) xtile(x, 5)))
Exp_hrd_scale=t(apply(Exp[signature_id,], 1, function(x) scale(x)))

####Test across race
wilcox.test(colMeans(Exp_hrd) ~ demo$Race, alternative='g')$p.value
wilcox.test(colMeans(Exp_hrd)[demo$Hist=='AD'] ~ demo$Race[demo$Hist=='AD'])$p.value
wilcox.test(colMeans(Exp_hrd)[demo$Hist=='SC'] ~ demo$Race[demo$Hist=='SC'])$p.value

wilcox.test(colMeans(Exp_hrd_scale) ~ demo$Race)$p.value
wilcox.test(colMeans(Exp_hrd_scale)[demo$Hist=='AD'] ~ demo$Race[demo$Hist=='AD'])$p.value
wilcox.test(colMeans(Exp_hrd_scale)[demo$Hist=='SC'] ~ demo$Race[demo$Hist=='SC'])$p.value
