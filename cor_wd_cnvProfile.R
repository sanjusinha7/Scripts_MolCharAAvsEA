##
library(ggpubr)
library(ggplot2)
cor_wd_cnvProfile<-function(vector=mat_CNV$chtp_Quan){
	corr_with_chtp_Quan = do.call(rbind,lapply(40:ncol(mat_CNV),function(x) 
				cor.test(vector, as.numeric(mat_CNV[,x]) )[c(3,4)]))
	rownames(corr_with_chtp_Quan)=colnames(mat_CNV)[40:ncol(mat_CNV)]
	corr_with_chtp_Quan=data.frame(corr_with_chtp_Quan)
	corr_with_chtp_Quan$fdr=p.adjust(corr_with_chtp_Quan[,1,], method='fdr')
	corr_with_chtp_Quan = corr_with_chtp_Quan[corr_with_chtp_Quan$fdr<0.1 & 
					corr_with_chtp_Quan$estimate>0.2,]
	corr_with_chtp_Quan=corr_with_chtp_Quan[order(unlist(corr_with_chtp_Quan$estimate), decreasing=T),]
	qq=data.frame(p.value=unlist(corr_with_chtp_Quan$p.value), 
		estimate=unlist(corr_with_chtp_Quan$estimate), fdr=unlist(corr_with_chtp_Quan$fdr))
	qq
}

cor_wd_CNV_burden=cor_wd_cnvProfile(mat_CNV$CNV_burden)
cor_wd_HRD_by_LOH=cor_wd_cnvProfile(mat_CNV$HRD_by_LOH)
cor_wd_chtp_Quan=cor_wd_cnvProfile(mat_CNV$chtp_Quan)


write.csv(qq, '/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/CNV_profile_corr_with_chtp_Quan.csv')



#########A few Genes
tiff('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/FBXW7_CNV.tiff')
theme_set(theme_bw(base_size = 20))
ggplot(mat_CNV, aes(y=FBXW7, x=race, fill=race))+
	geom_boxplot()+
	geom_jitter(position=position_jitter(0.05))+
	stat_compare_means(label = "p.signif")+
        facet_grid( ~hist, scales = "free", space = "free_x")+
	ggtitle('FBXW7 CNV Profile across Races')

dev.off()
mat_CNV_xtile=mat_CNV

posNeg<-function(x){
	if(x>0) {'Gain'
	} else if(x<0) {'loss'
	} else if(x==0) {'Neutral' }
}

mat_CNV_xtile[,-(1:39)] = apply(mat_CNV[,-(1:39)], 2, function(x) sapply(x, posNeg)  )

printPlot<-function(geneName, fit1){
	tiff(paste('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/',geneName,'_KM.tiff', sep=''))
	theme_set(theme_bw(base_size = 20))
	print(ggsurvplot_facet(fit1, mat_CNV_xtile, facet.by = "hist", pval=TRUE, pallette='jco'))
	dev.off()
}
GeneList=c('BOP1','DGAT1','HSF1','CYP2E1','FBXW7')

fit=list(fit1=survfit( Surv(survival, lungcancer_death_all_years) ~ BOP1, mat_CNV_xtile ),
fit2=survfit( Surv(survival, lungcancer_death_all_years) ~ DGAT1, mat_CNV_xtile ),
fit3=survfit( Surv(survival, lungcancer_death_all_years) ~ HSF1, mat_CNV_xtile ),
fit4=survfit( Surv(survival, lungcancer_death_all_years) ~ CYP2E1, mat_CNV_xtile ),
fit5=survfit( Surv(survival, lungcancer_death_all_years) ~ FBXW7, mat_CNV_xtile ))


lapply(1:length(fit), function(x) printPlot(GeneList[x], fit[[x]]))


