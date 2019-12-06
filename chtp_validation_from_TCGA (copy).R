##required libraries
require(parallel)
require(ggpubr)
require(ggplot2)

###Validation in TCGA::Higher Chromothripsis in LUSC-AA
#load('/cbcb/project2-scratch/sanju/ISLE1.0/TCGA.RData')
load('/Users/sinhas8/Downloads/TCGA.RData')
prob$race=as.factor(prob$race)
levels(prob$race)[c(3,6)]=c('AA', 'EA')

#source('/home/sinhas8/Projects/Project_Chromotrypsis/3.Tools/BPC_Test.R')
setwd('/Users/sinhas8/Project_Chromotrypsis/2.Data/Seg_TCGA')
Cancer_Types_list=list.files()
source('/Users/sinhas8/Project_Chromotrypsis/3.Tools/BPC_Test.R')

############################################################################################################################################################################################
############################################################################################################################################################################################
find_chtp<-function(Cancer_Type='LUSC'){
	if(sum(grepl(Cancer_Type, prob$types))>0){
		seg=read.csv(paste(Cancer_Type,list.files(Cancer_Type), sep='/'), sep='\t') 
		##Preprocessing
		colnames(seg)=c('File','chr','start', 'end', 'no_of_markers','CNV')
		seg_wdout_Diplod=seg[seg$CNV < -0.1 | seg$CNV >0.1 ,]
		seg_wdout_Diploid_List=split(seg_wdout_Diplod, seg_wdout_Diplod$File)
		
		Sample_Type=sapply(names(seg_wdout_Diploid_List), function(x) 
                                   paste(strsplit(x, '-')[[1]][4], collapse='-') )
		Tumor_samples = names(seg_wdout_Diploid_List)[which(as.numeric(substring(Sample_Type, 1,2))<10)]
	
		seg_wdout_Diploid_List=seg_wdout_Diploid_List[Tumor_samples]
		names(seg_wdout_Diploid_List)  = sapply(Tumor_samples, function(x) 
                                                 paste(strsplit(x, '-')[[1]][1:3], collapse='-') )
		##Testing Hypothesis 1
		after_H1=sapply(seg_wdout_Diploid_List, function(x) 
                          split(x, x$chr)[sapply(split(x, x$chr), nrow)>mean(sapply(split(x, x$chr), nrow))] )

		##H2: Iterative cluster	##Functions needed
		qq = lapply(after_H1, function(x) BPC_Test(sample=x))
		qq=qq[!sapply(qq, function(x) is.null(nrow(x)))]

		Significance     = sapply(qq, function(x) unlist(x[1,]))
		Number_of_Events = sapply(qq, function(x) unlist(x[3,]))
		Median_Distance  = sapply(qq, function(x) unlist(x[4,]))

		Samples_wd_Chtp=sapply( Significance, function(x) names(which(p.adjust(x, method='fdr')<0.1)) )
		Samples_wd_Chtp=Samples_wd_Chtp[sapply(Samples_wd_Chtp, length)>0]
		Samples_wd_Chtp=data.frame(Sample_names=as.character(names(Samples_wd_Chtp)),
                                           chromosome=paste(Samples_wd_Chtp))
		Samples_wd_Chtp[,1]= as.character(Samples_wd_Chtp[,1])

		#Quantity
		chtp_Quantification=sapply(Samples_wd_Chtp$chromosome, function(x) 
                                     length(strsplit(gsub('\\"','',gsub('[c\\(\\)]','',x)), ',')[[1]]) )
		Samples_wd_Chtp$Quantity=chtp_Quantification

		#Post-hoc analysis		#histogram for later
#		print(sort(table(unlist(sapply(gsub('[\\"\\(\\)c]','', Samples_wd_Chtp$chromosome), function(x) strsplit(x, ', '))))))

		chtp_df=data.frame(sample_name=names(seg_wdout_Diploid_List),
                                   chtp=as.numeric(!is.na(match(names(seg_wdout_Diploid_List),
                                   Samples_wd_Chtp$Sample_names))))
		chtp_df$Chr=''
		chtp_df$Chr[!is.na(match(names(seg_wdout_Diploid_List),
		                         Samples_wd_Chtp$Sample_names))] =as.character(Samples_wd_Chtp$chromosome)
		chtp_df_matched=chtp_df[match(prob$samples[prob$types==Cancer_Type], chtp_df$sample_name),]

		chtp_df_matched=data.frame(chtp_df_matched, survival= prob$surv.dt$time[prob$type==Cancer_Type],
	                      lungcancer_death_all_years=prob$surv.dt$status[prob$type==Cancer_Type],
	                      race=prob$race[prob$type==Cancer_Type], hist=Cancer_Type )
	
		chtp_df_matched=chtp_df_matched[which(chtp_df_matched$race=='AA'|
                                                      chtp_df_matched$race=='EA'),]
	} else(chtp_df_matched=NA)
		chtp_df_matched
}

##
CNV_Types_Limit=1
TCGA_chtp=mclapply(1:length(Cancer_Types_list), function(x) err_handle(find_chtp(Cancer_Type=Cancer_Types_list[x])), 
                                                            mc.cores=detectCores())
saveRDS(TCGA_chtp, '/cbcb/project2-scratch/sanju/Chromotrypsis/2.Data/TCGA_chtp.RDS')
#TCGA_chtp = readRDS('/cbcb/project2-scratch/sanju/Chromotrypsis/2.Data/TCGA_chtp.RDS')
#saveRDS(TCGA_chtp, '/cbcb/project2-scratch/sanju/Chromotrypsis/2.Data/TCGA_chtp_UpdatedDef.RDS')

TCGA_chtp=readRDS('/Users/sinhas8/Project_Chromotrypsis/2.Data/TCGA_chtp_UpdatedDef.RDS')
head(TCGA_chtp)
TCGA_chtp_FF=TCGA_chtp
TCGA_chtp_FF=TCGA_chtp[!sapply(TCGA_chtp, function(x) is.null(dim(x)))]
TCGA_chtp_df=TCGA_chtp_FF
TCGA_chtp_df=do.call(rbind, TCGA_chtp_FF)
TCGA_chtp_df_filtered=TCGA_chtp_df[TCGA_chtp_df$race=='AA'|TCGA_chtp_df$race=='EA',]
TCGA_chtp_df_filtered$race=factor(TCGA_chtp_df_filtered$race)
TCGA_chtp_df_filtered=na.omit(TCGA_chtp_df_filtered)

head(TCGA_chtp_df_filtered)
##
write.csv(TCGA_chtp_df_filtered, '/Users/sinhas8/Project_Chromotrypsis/2.Data/TCGA_chtp_with_CHR.csv', 
          quote=F, row.names=F)

# Start from here
TCGA_chtp_df_filtered=read.csv('/Users/sinhas8/Project_Chromotrypsis/2.Data/TCGA_chtp_with_CHR.csv')

AA_Count_greater_atleast5 = sapply(split(TCGA_chtp_df_filtered$race, 
                                         TCGA_chtp_df_filtered$hist), 
                                   function(x) table(x)[table(x)>0][1]>4)
table(TCGA_chtp_df_filtered$hist)
TCGA_chtp_df_filtered =do.call(rbind, split(TCGA_chtp_df_filtered,
                                            TCGA_chtp_df_filtered$hist)[AA_Count_greater_atleast5])
df=TCGA_chtp_df_filtered

###
chtp_distribution<-function(cancer_mat=temp[[2]]){
	chtp_dist=data.frame(table(unlist(sapply(cancer_mat$Chr, function(x) 
                     strsplit(gsub('[c\\(\\) \\"]','',x), '_')[[1]]) ))/nrow(cancer_mat))
	chtp_dist$Var1=sapply(chtp_dist$Var1, function(x) as.numeric(as.character(x)))
	chtp_dist=chtp_dist[chtp_dist$Var1!='24',]
	if(nrow(chtp_dist)<23){
	Rest=data.frame(Var1=(1:23)[-na.omit(match(chtp_dist$Var1, 1:23))], Freq=0)
	data.frame(rbind(Rest, chtp_dist)[order(as.numeric(rbind(Rest, chtp_dist)[,1]), decreasing=F),], 
	           Cancer_Type=cancer_mat$hist[1], 
	           Race=cancer_mat$race[1])
	} else {
	  data.frame(chtp_dist[order(as.numeric(chtp_dist[,1]), decreasing=F),], 
	             Cancer_Type=cancer_mat$hist[1],
	             Race=cancer_mat$race[1])
	}
	
}
df$hist=as.character(df$hist)
temp=split(df, list(df$hist,df$race))
df_list=lapply(1:length(temp), function(x) err_handle(chtp_distribution(temp[[x]])) )
df2plot=do.call(rbind, df_list)
ID <- 1:23
df2plot$Var1=as.numeric(df2plot$Var1)

levels(df2plot$Cancer_Type)
colorType_Set='Set1'
# tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/CHTP_9thMay', height=2000, width=3000)
Ext_Fig2C=ggplot(df2plot, aes(y=Freq, x=Var1, fill=Race))+
  geom_bar(stat = "identity", position=position_dodge())+
  facet_wrap(~Cancer_Type, scales='free')+
  labs(y='Chromothripsis Frequency', x='Chromosomes')+
  scale_x_continuous('Chromosomes', labels = as.character(ID), breaks = ID)+
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(size = 8))+
  scale_fill_brewer(palette=colorType_Set)+
  guides(fill=FALSE)
# dev.off()

#########

# ad_AA = which(mat_filtered$hist=='adeno' & mat_filtered$race=='AA' )
# ad_EA = which(mat_filtered$hist=='adeno' & mat_filtered$race=='EA' )
# sq_AA = which(mat_filtered$hist=='sq' & mat_filtered$race=='AA' )
# sq_EA = which(mat_filtered$hist=='sq' & mat_filtered$race=='EA' )


mat$chtp
chtp_distribution<-function(Type="LUSC", test_race='AA'){
  chtp_dist=data.frame(table(unlist(sapply(mat$chtp[mat$hist==Type & mat$race==test_race], function(x) 
    strsplit(gsub('[c\\(\\) \\"]','',x), ',')[[1]]) ))/length(Type))
  Rest=data.frame(Var1=(1:24)[-match(chtp_dist$Var1, 1:24)], Freq=0)
  rbind(Rest, chtp_dist)[order(as.numeric(rbind(Rest, chtp_dist)[,1]), decreasing=F),]
}

df=do.call(rbind,list(data.frame(chtp_distribution(Type="LUAD", test_race='AA'), Hist='LUAD', Race='AA'),
                      data.frame(chtp_distribution(Type="LUAD", test_race='EA'), Hist='LUAD', Race='EA'),
                      data.frame(chtp_distribution(Type="LUSC", test_race='AA'), Hist='LUSC', Race='AA'),
                      data.frame(chtp_distribution(Type="LUSC", test_race='EA'), Hist='LUSC', Race='EA')))
df$Var1=sapply(df$Var1, as.numeric)
df=df[!(df$Var1== '23' | df$Var1== '24'),]

ID <- 1:22


Ext_Fig2B=ggplot(df, aes(y=Freq, x=Var1, fill=Race))+
  geom_bar(stat = "identity", position=position_dodge())+
  facet_wrap(~Hist, scales = 'free', shrink = T, nrow=2)+
  labs(x='Chromosomes', y='Chromothripsis Frequency')+
  scale_x_continuous("ID", labels = as.character(ID), breaks = ID)+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(size = 10))+
  scale_fill_brewer(palette=colorType_Set)


tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Extended_Figure2.tif', height=1800, width=2200)
grid.arrange(plot_grid(Ext_Fig2A, labels = 'A', label_size = 30),
             plot_grid(Ext_Fig2B, labels = 'B', label_size = 30),
             plot_grid(Ext_Fig2C, labels = 'C', label_size = 30),
          layout_matrix = rbind(c(1, 1, 2),
                                c(3, 3, 3),
                                c(3, 3, 3),
                                c(3, 3, 3))
)
dev.off()
