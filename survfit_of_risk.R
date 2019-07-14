require(survival)

qq=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/mat_2ndOct.csv', sep='\t')
qq$CNV_burden=GI
qq$extreme_cnv= GI>0.75 | GI<0.25
write.csv(qq,'/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/mat_18thOct.csv')


fit1<-survfit(Surv(survival, lungcancer_death_all_years) ~ extreme_cnv+race, qq[sq,])
fit2<-survfit(Surv(survival, lungcancer_death_all_years) ~ extreme_cnv+race, qq[ad,])

p1<-ggsurvplot(fit1, data = qq[sq,], pval=TRUE, title='For NSCLC subtype Squamous Cell Carcinoma', legend='none', risk.table = TRUE)
p2<-ggsurvplot(fit2, data = qq[ad,], pval=TRUE, title='For NSCLC subtype Adenocarcinoma', legend='right', risk.table = TRUE)

tiff('/home/sinhas8/Projects/Project_Chromotrypsis/Results_New/Kaplan_Mier_GI_extrem.tiff', height=800, width=1600)
arrange_ggsurvplots(list(p1, p2), ncol=2, nrow=1)
dev.off()


###Cox regression
coxph(Surv(survival, lungcancer_death_all_years) ~ extreme_cnv+race+stage, qq[sq,])
coxph(Surv(survival, lungcancer_death_all_years) ~ extreme_cnv+race+stage, qq[ad,])




