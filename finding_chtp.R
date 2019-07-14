###Finding Chromothripsis::
require(ggplot2)
require(ggpubr)
require(parallel)
#Preprocessing
seg=read.table('/Users/sinhas8/Downloads/Both_Phases_combined.txt')
colnames(seg)=c('File','chr','start', 'end', 'no_of_markers','CNV')
#removing diploid-regions
seg_wdout_Diploid=seg[seg$CNV!=2,]

seg_wdout_Diploid_List=split(seg_wdout_Diploid, seg_wdout_Diploid$File)

############################################################################################################################################################################################
#Hypothesis 1:: Signal-background Hypothesis::The alteration events on a Chromotripsis Chromosome(CC) should be higher than the background events. 
############################################################################################################################################################################################
after_H1=sapply(seg_wdout_Diploid_List, function(x) split(x, x$chr)[sapply(split(x, x$chr), nrow)>mean(sapply(split(x, x$chr), nrow))] )

##H2: Iterative cluster
##Functions needed
source('/Users/sinhas8/Project_Chromotrypsis/3.Tools/BPC_Test.R')
CNV_Types_Limit=1
#mat=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/mat_2ndOct.csv', sep='\t')
qq=mclapply(after_H1, function(x) BPC_Test(sample=x), mc.cores=detectCores())
qq=qq[!sapply(qq, function(x) is.null(nrow(x)))]

Significance     = sapply(qq, function(x) unlist(x[1,]))
Number_of_Events = sapply(qq, function(x) unlist(x[3,]))
Median_Distance  = sapply(qq, function(x) unlist(x[4,]))

Samples_wd_Chtp=sapply( Significance, function(x) names(which(p.adjust(x, method='fdr')<0.1)) )
Samples_wd_Chtp=Samples_wd_Chtp[sapply(Samples_wd_Chtp, length)>0]
Samples_wd_Chtp=data.frame(Sample_names=as.character(names(Samples_wd_Chtp)), chromosome=paste(Samples_wd_Chtp))
Samples_wd_Chtp[,1]= as.character(Samples_wd_Chtp[,1])

#Quantity
chtp_Quantification=sapply(Samples_wd_Chtp$chromosome, function(x) length(strsplit(gsub('\\"','',gsub('[c\\(\\)]','',x)), ', ')[[1]]) )
Samples_wd_Chtp$Quantity=chtp_Quantification


##adding inconsistent names::
Inconsistent_Names=Samples_wd_Chtp$Sample_names[is.na(match(sapply(Samples_wd_Chtp$Sample_names,function(x) gsub('-','\\.', x)), sapply(mat$File,function(x) gsub('-','\\.', x))))]
Correct_Names = sapply(sapply(Inconsistent_Names, function(x) strsplit(strsplit(as.character(x), '\\.')[[1]][1], '_')[[1]][2]), function(x) grep(x, mat$names, value=T))


Samples_wd_Chtp$Sample_names[is.na(match(sapply(Samples_wd_Chtp$Sample_names,function(x) gsub('-','\\.', x)), sapply(mat$File,function(x) gsub('-','\\.', x))))]= Correct_Names

##Part 2 due to inconsistent part##
mat$chtp_Quan_new=0; mat$chtp_new=''
mat$chtp_Quan_new[match(Samples_wd_Chtp$Sample_names, mat$File)] = Samples_wd_Chtp$Quantity
mat$chtp_new[match(Samples_wd_Chtp$Sample_names, mat$File)] = as.character(Samples_wd_Chtp$chromosome)

#mat_filtered=mat[mat$hist=='sq' | mat$hist=='adeno',]

mat_filtered$chromothripsis_presence_new= as.numeric(as.logical(mat_filtered$chtp_Quan_new))

means <- aggregate(chromothripsis_presence_new ~  race+hist, mat_filtered, mean)
print(means)

tiff('/home/sinhas8/Projects/Project_Chromotrypsis/Results_New/Chtp_MeanChromosomes_violen_delthis_New_wd1Type.tiff')
ggplot(mat_filtered, aes(y=chromothripsis_presence_new, x=race, fill=race))+
	geom_violin()+
#	coord_flip()+
#	geom_boxplot(width=0.1)+
	facet_wrap(~hist)+
	stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
	stat_compare_means(label = "p.signif", method = "wilcox.test",  label.y = 1.2,  label.x = 1.5)+
	ylab('Mean Number of Chromosomes')

#	  theme(plot.background = element_blank(), 
#	   panel.grid.major = element_blank(),
#	   panel.grid.minor = element_blank(), 
#	   panel.border = element_blank(),
#	   panel.background = element_blank(),
#	   axis.title.x = element_blank(),
#	   axis.title.y = element_blank(),
#	   axis.text.x = element_blank(), 
#	   axis.text.y = element_blank(),
#	   axis.ticks = element_blank()
#	     )
dev.off()

##No significant effect
#summary(coxph(Surv(mat$survival, mat$lungcancer_death_all_years)~ as.logical(mat$chtp_Quan))))

#Chr 8 effect
#mat$chr_eight_chtp=NA
#mat$chr_eight_chtp[as.logical(mat$chtp_Quan)]=0
#mat$chr_eight_chtp[which(mat$chtp=='8' | grepl('*8,*',mat$chtp))]=1
ad_AA = which(mat_filtered$hist=='adeno' & mat_filtered$race=='AA' )
ad_EA = which(mat_filtered$hist=='adeno' & mat_filtered$race=='EA' )
sq_AA = which(mat_filtered$hist=='sq' & mat_filtered$race=='AA' )
sq_EA = which(mat_filtered$hist=='sq' & mat_filtered$race=='EA' )

chtp_distribution<-function(Type){
	chtp_dist=data.frame(table(unlist( sapply(mat_filtered$chtp_new[Type], function(x) 
                             strsplit(gsub('[c\\(\\) \\"]','',x), ',')[[1]]) ))/length(Type))
	Rest=data.frame(Var1=(1:24)[-match(chtp_dist$Var1, 1:24)], Freq=0)
	rbind(Rest, chtp_dist)[order(as.numeric(rbind(Rest, chtp_dist)[,1]), decreasing=F),]
}

df=do.call(rbind,list(data.frame(chtp_distribution(ad_AA), Hist='ADC', Race='AA'),
data.frame(chtp_distribution(ad_EA), Hist='ADC', Race='EA'),
data.frame(chtp_distribution(sq_AA), Hist='SCC', Race='AA'),
data.frame(chtp_distribution(sq_EA), Hist='SCC', Race='EA')))
df$Var1=sapply(df$Var1, as.numeric)

ID <- 1:24

tiff('/home/sinhas8/Projects/Project_Chromotrypsis/Results_New/Chtp_Chr_dist_our_Cohort.tiff', height=300, width=600)
theme_set(theme_bw(base_size = 20))
ggplot(df, aes(y=Freq, x=Var1, fill=Race))+
	geom_bar(stat = "identity", position=position_dodge())+
#	coord_flip()+
#	geom_boxplot(width=0.1)+
	facet_wrap(~Hist, nrow=2, ncol=1)+
#	stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
#              width = 0.75, size = 1, linetype = "solid")+
#	stat_compare_means(label = "p.signif", method = "wilcox.test",  label.y = 1.2,  label.x = 1.5)+
	labs(y='Chromothripsis Frequency', x='Chromosomes')+
	scale_x_continuous("ID", labels = as.character(ID), breaks = ID)+
	theme(axis.text.x = element_text(size = 10))
#	  theme(plot.background = element_blank(), 
#	   panel.grid.major = element_blank(),
#	   panel.grid.minor = element_blank(), 
#	   panel.border = element_blank(),
#	   panel.background = element_blank(),
#	   axis.title.x = element_blank(),
#	   axis.title.y = element_blank(),
#	   axis.text.x = element_blank(), 
#	   axis.text.y = element_blank(),
#	   axis.ticks = element_blank()
#	     )
dev.off()

