###validation in tcga::higher genomic instability
#libraies needed
require(survival)
require(survminer)
require('doMC')
require(ggplot2)
registerDoMC(cores=detectCores())

##files needed
setwd('/cbcb/project2-scratch/sanju')
Cancer_Types=list.files()

##loading tcga data
load('/cbcb/project2-scratch/sanju/ISLE1.0/TCGA.RData')
prob$race=as.factor(prob$race)
levels(prob$race)[3]='AA'
levels(prob$race)[6]='EA'

##
processing2gi<-function(cancer_type, event_type='GI'){
  type=cancer_type
  if(sum(grepl(type, prob$types))>0){
		print(type)
		mat=read.csv(paste(type,list.files(type), sep='/'), sep='\t')
		##preprocessing:: 
		colnames(mat)=c('sample', 'chromosome', 'start', 'end', 'num_probes', 'segment_mean')
	  if(event_type=='GI'){
		mat_wdout_diplod=mat[mat$segment_mean < -0.1 | mat$segment_mean >0.1 ,]
		} else if(event_type=='loss_GI'){
		mat_wdout_diplod=mat[mat$segment_mean < -0.1,]
		} else if(event_type=='gain_GI'){
		mat_wdout_diplod=mat[mat$segment_mean > 0.1,]
		}
		##considering non-diploid regions
		mat_list=split(mat_wdout_diplod, mat_wdout_diplod$sample)
	
		gi_tcga=lapply(mat_list, function(x) sum(unlist(lapply(split(x, x$chromosome), 
	                       function(y) apply(y, 1, function(z) as.numeric(z[4]) -as.numeric(z[3]) ) ))  ) )
		gi_tcga=unlist(gi_tcga)/3.3e9
	
		sample_type=sapply(names(gi_tcga), function(x) paste(strsplit(x, '-')[[1]][4], collapse='-') )
		tumor_samples=which(as.numeric(substring(sample_type, 1,2))<10)

		gi_tumor=gi_tcga[tumor_samples]
		names(gi_tumor)  = sapply(names(gi_tumor), function(x) 
                                          paste(strsplit(x, '-')[[1]][1:3], collapse='-') )
		gi_tumor_matched = gi_tumor[match(prob$samples[ prob$type==type], names(gi_tumor))]
		extreme_gi=gi_tumor_matched>0.75 | gi_tumor_matched<0.25 
		qq=data.frame(gi=gi_tumor_matched, extreme_cnv=extreme_gi, 
                              survival= prob$surv.dt$time[prob$type==type],
	                      lungcancer_death_all_years=prob$surv.dt$status[prob$type==type],
	                      race=prob$race[prob$type==type], hist=type )		
#		obj_return=qq
		obj_return=qq[which(qq$race=='AA'|qq$race=='EA'),]
	} else{obj_return=NA}
	obj_return
}

##
qq<-foreach(i= seq(Cancer_Types)) %dopar% {
	tryCatch(processing2gi(cancer_type=Cancer_Types[i],
	                       event_type='loss_GI'), error=function(err){NULL})
}
names(qq)=Cancer_Types

#Preprocessing
qq=qq[names(unlist(sapply(qq, nrow))>0)]
qq_df=do.call(rbind, qq)
qq_df=na.omit(qq_df)

##cancer types with at least 5 samples
Race_Distribution = sapply(split(qq_df$race, qq_df$hist), function(x) table(x)[table(x)>0])
AA_Count_greater_atleast5 = sapply(split(qq_df$race, qq_df$hist), function(x) table(x)[table(x)>0][1]>4)

#qq_df_sorted=do.call(rbind, split(qq_df, qq_df$hist)[GI_Median_Order])
qq_df =do.call(rbind, split(qq_df, qq_df$hist)[AA_Count_greater_atleast5])

##Median order
GI_Median=sapply(split(qq_df$gi, qq_df$hist), median)
GI_Median_Order=order(GI_Median)
qq_df$hist<-factor(qq_df$hist, levels=levels(qq_df$hist)[GI_Median_Order])
GI_Median=data.frame(GI_Median)

##Tissue_Type
qq_df$Tissue_Type = 'Ungrouped'
##All hemotalogical and lymphatic 
qq_df$Tissue_Type[  qq_df$hist=='THYM' |  qq_df$hist=='DLBC'|  qq_df$hist=='LAML'] = 'Hema&lymph'
##All Solid
qq_df$Tissue_Type[  qq_df$hist=='OV' |  qq_df$hist=='UCEC'  | qq_df$hist=='CESC' |  qq_df$hist=='BRCA' |
                    qq_df$hist=='BLCA' | qq_df$hist=='PRAD' | qq_df$hist=='TGCT' |
                    qq_df$hist=='KIRC' | qq_df$hist=='KICH' | qq_df$hist=='KIRP' | 
		    qq_df$hist=='THCA' | qq_df$hist=='ACC'  |
                    qq_df$hist=='ESCA' | qq_df$hist=='STAD' | qq_df$hist=='COAD' | qq_df$hist=='READ'  |
                    qq_df$hist=='LIHC' | qq_df$hist=='PAAD' | qq_df$hist=='CHOL' |
		    qq_df$hist=='HNSC' |
                    qq_df$hist=='LUAD' | qq_df$hist=='LUSC' | qq_df$hist=='MESO'] = 'Solid'
#Neural-crest-derived
qq_df$Tissue_Type[  qq_df$hist=='GBM'  | qq_df$hist=='LGG'  |
                    qq_df$hist=='SARC' | qq_df$hist=='UCS'  | 
		    qq_df$hist=='PCPG' | qq_df$hist=='SKCM' |
		    qq_df$hist=='UVM'] = 'Neural-crest-derived'

##CEll_Type
qq_df$Cell_Type = 'Ungrouped'
qq_df$cell_Type[  qq_df$hist=='OV' |  qq_df$hist=='UCEC'  | qq_df$hist=='CESC' |  qq_df$hist=='BRCA']='gynecologic'
qq_df$cell_Type[qq_df$hist=='BLCA' | qq_df$hist=='PRAD' | qq_df$hist=='TGCT' |
                qq_df$hist=='KIRC' | qq_df$hist=='KICH' | qq_df$hist=='KIRP']='Urological'
qq_df$cell_Type[qq_df$hist=='THCA' | qq_df$hist=='ACC']='Endocrine'
qq_df$cell_Type[qq_df$hist=='ESCA' | qq_df$hist=='STAD' | qq_df$hist=='COAD' | qq_df$hist=='READ']='core-GI'
qq_df$cell_Type[qq_df$hist=='LIHC' | qq_df$hist=='PAAD' | qq_df$hist=='CHOL']='Dev-GI'
qq_df$cell_Type[qq_df$hist=='HNSC' | qq_df$hist=='LUAD' | qq_df$hist=='LUSC' | qq_df$hist=='MESO']='thoracic-H&N'

qq_df$cell_Type[qq_df$hist=='GBM'  | qq_df$hist=='LGG']='CNS'
qq_df$cell_Type[qq_df$hist=='SARC' | qq_df$hist=='UCS' | qq_df$hist=='PCPG']='Soft' 
qq_df$cell_Type[qq_df$hist=='SKCM' | qq_df$hist=='UVM']='Skin-Eye' 
##
qq_df$CellofOrigin='Rest'
qq_df$CellofOrigin[ qq_df$hist=='LUSC' |  qq_df$hist=='ESCA' |  qq_df$hist=='HNSC'|  qq_df$hist=='CESC'|  qq_df$hist=='BLCA']='Pan-Squamous'
qq_df$CellofOrigin[ qq_df$hist=='PRAD' |  qq_df$hist=='STAD' |  qq_df$hist=='COAD'|  qq_df$hist=='READ'|  qq_df$hist=='PAAD'|  qq_df$hist=='LUAD']='Pan-Adeno'
qq_df$CellofOrigin[ qq_df$hist=='KIRC' |  qq_df$hist=='KIRP']='Pan-Kidney'

qq_df_gain=qq_df
qq_df_loss=qq_df
df=data.frame(qq_df_gain, GI_gain=qq_df_loss[,1])
write.csv(df, '/home/')

#wilcox.test(qq_df$gi ~ qq_df$race, alternative='g')
##
tiff('/cbcb/project2-scratch/sanju/tcga_boxplot_byrace_cnv_6thDec_burden_acrossraces_inhist_19th_atleast4.tiff' , height=1000, width=2000)
theme_set(theme_bw(base_size = 20))
ggplot(qq_df, aes(y=gi, x=hist, color=race))  +
	geom_boxplot()+
	stat_compare_means(method = 'wilcox', label = "p.signif")+
        facet_grid( ~Tissue_Type + CellofOrigin, scales = "free", space = "free_x")
dev.off()

#tiff('/cbcb/project2-scratch/sanju/delthis_Gain.tiff' , height=500, width=1000)
#plot_grid(GI_plot, Race_dist_plot, ncol=1, nrow=2, rel_heights =c(9/10,1/10), align = 'h')
#dev.off()

##
tiff('/cbcb/project2-scratch/sanju/tcga_dist_of_cnvGain_burden_acrossraces_inhist_19th_New.tiff', height=1200, width=1000)
theme_set(theme_bw(base_size = 20))
ggplot(qq_df, aes(gi, color=race))  +
  geom_density() +
  geom_vline(xintercept=0.25, linetype="dashed", color = "black")+  
  geom_vline(xintercept=0.75, linetype="dashed", color = "black")+
#  facet_grid(~ , scales = "free", space = "free_x")
  facet_wrap(CellofOrigin ~ hist,  nrow = 6, ncol = 5, scales = "free")

#  geom_vline(xintercept=mean(gi), linetype="dashed", color = "black")+  
# geom_vline(xintercept=0.75, linetype="dashed", color = "black")+

dev.off()

###survival 
p<-foreach(i= seq(length(qq))) %dopar% {
	fit<-survfit(Surv(survival, lungcancer_death_all_years) ~ extreme_cnv, qq[[i]])
	ggsurvplot(fit, data = qq[[i]], pval=TRUE,
                   title=names(qq[i]), 
#		   legend='none',
		   risk.table = FALSE)

}

##
tiff('/cbcb/project2-scratch/sanju/kaplan_mier_tcga_filtered_only_GI.tiff', height=1250, width=1500)
arrange_ggsurvplots(p, ncol=6, nrow=5)
dev.off()
pan_sq
df2plot=data.frame(Sig=p.adjust(tt1$Sig, method='fdr'), AA_proportion=tt1$AA_Proportion)
df2plot$whether_Significant=0
df2plot$whether_Significant[df2plot$Sig<0.1]=1
df2plot$whether_Significant= factor(df2plot$whether_Significant, labels=c('Not Significant', 'Significant'))

#P-value vs proportion of GI
tiff('/cbcb/project2-scratch/sanju/Chromotrypsis/Pvalue_vs_AAProportion_PanSq.tiff')
theme_set(theme_bw(base_size = 20))
ggplot(df2plot[pan_sq,], aes(x=AA_proportion, y=Sig, color=whether_Significant))+
	geom_point()+
#	annotate("text", x=0.35, y=0.8, label= "Spearman\n Rho= -0.35,\n P<0.14",  size = 5)+ 
	annotate("text", x=0.1, y=0.25, label= "Spearman Rho= -0.66, P<0.21\nPearson rho= -0.66, P<0.02 ", size = 5)+ 
#	facet_wrap(~Hist, nrow=2, ncol=1)+
#	stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
#              width = 0.75, size = 1, linetype = "solid")+
#	stat_compare_means(label = "p.signif", method = "wilcox.test",  label.y = 1.2,  label.x = 1.5)+
	labs(y='FDR Corrected', x='AA Samples Proportion')+
	theme(legend.title=element_blank(), legend.position="top")
#	theme(axis.text.x = element_text(size = 10))
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

