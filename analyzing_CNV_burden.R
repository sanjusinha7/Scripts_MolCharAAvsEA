##In this script our aim to reproduce the results of Mariam Jamal-Hanjani et al of NJC where they have  showed subclonal CNV proportions is a prognostic marker and is associated with ITH in adenocarcinoma.


demo=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/demo_and_clones.csv')
demo$stage_trimmed=demo$stage

levels(demo$stage_trimmed)=c('NA','1','1','1','2','3','1','1','2','2','3','3')

clone_p1=read.csv('/home/sinhas8/Downloads/P2_Initial_Data/OncoClone_CloneInfo.csv')
clone_p1=clone_p1[clone_p1$File != '',]
clone_p2=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/OncoClone_Metric_Output.TXT', sep='\t')
clone_p2=clone_p2[clone_p2$File != '',]
colnames(clone_p1)=colnames(clone_p2)
clone_info=rbind(clone_p1, clone_p2)

demo$names=as.character(demo$names)
clone_info=clone_info[match(as.character(demo$names), c(as.character(clone_p1$File), as.character(clone_p2$File))),]

mat=cbind(demo, clone_info)
mat$purity=as.numeric(gsub('%','',mat$X.AC))
mat$CNV_burden=new_CNV_burden
write.csv(mat, '/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/demo_and_clones.csv')



##Whether the parameters ahs any clinical relevance

clinical_relevance<-function(factor=mat$Clones, histtype='adeno'){
	summary(coxph(Surv(mat$survival[mat$hist==histtype], mat$lungcancer_death_all_years[mat$hist==histtype]) ~ factor[mat$hist==histtype] + mat$stage_trimmed[ mat$hist==histtype]+ mat$race[ mat$hist==histtype] ))$coefficients[c(1,5),c(1,5)]
}


###Different info we have about sample::
# Histology/ Race/ Stage/ Survival/ :: Purity/ Ploidy/ Static CNV_Burden/  Clonality::

##Purity vs histology race
wilcox.test(mat$purity[mat$hist=='adeno'] ~ mat$race[mat$hist=='adeno'])
wilcox.test(mat$purity[mat$hist=='sq'] ~ mat$race[mat$hist=='sq'])

new_CNV_burden=readRDS('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/delthis.RDS')
