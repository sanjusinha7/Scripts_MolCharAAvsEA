##
df_adeno$Del_Sig='Recurrent in none'
df_adeno$Del_Sig[df_adeno$AA_Del_qvalue<0.1 & df_adeno$EA_Del_qvalue<0.1]='Recurrent in Both'
df_adeno$Del_Sig[df_adeno$AA_Del_qvalue<0.1 & df_adeno$EA_Del_qvalue>0.1]='Recurrent in AA Only'
df_adeno$Del_Sig[df_adeno$AA_Del_qvalue>0.1 & df_adeno$EA_Del_qvalue<0.1] = 'Recurrent in EA Only'
table(df_adeno$Del_Sig)

df_adeno$Amp_Sig='Recurrent in none'
df_adeno$Amp_Sig[df_adeno$AA_Amp_qvalue<0.1 & df_adeno$EA_Amp_qvalue<0.1]='Recurrent in Both'
df_adeno$Amp_Sig[df_adeno$AA_Amp_qvalue<0.1 & df_adeno$EA_Amp_qvalue>0.1]='Recurrent in AA Only'
df_adeno$Amp_Sig[df_adeno$AA_Amp_qvalue>0.1 & df_adeno$EA_Amp_qvalue<0.1] = 'Recurrent in EA Only'
table(df_adeno$Amp_Sig)

##
df1_Amp=df1_Amp[!is.na(df1_Amp$Exp_Assoc_Sig),]
df1_Amp_Sig_sq=df1_Amp[p.adjust(df1_Amp$Exp_Assoc_Sig, method='fdr')<0.1,]
df1_Del=df1_Del[!is.na(df1_Del$Exp_Assoc_Sig),]
df1_Del_Sig_sq=df1_Del[p.adjust(df1_Del$Exp_Assoc_Sig, method='fdr')<0.1,]

##
df1_Amp=df1_Amp[!is.na(df1_Amp$Exp_Assoc_Sig),]
df1_Amp_Sig_ad=df1_Amp[df1_Amp$Exp_Assoc_Sig<0.1,]
df1_Del=df1_Del[!is.na(df1_Del$Exp_Assoc_Sig),]
df1_Del_Sig_ad=df1_Del[df1_Del$Exp_Assoc_Sig<0.1,]


##
cat(c(rownames(df1_Del_Sig_ad[df1_Del_Sig_ad$Aberrance_EA=='',]), rownames(df1_Amp_Sig_ad[df1_Amp_Sig_ad$Aberrance_EA=='',])), sep='\n')
cat(c(rownames(df1_Del_Sig_sq[df1_Del_Sig_sq$Aberrance_EA=='',]), rownames(df1_Amp_Sig_sq[df1_Amp_Sig_sq$Aberrance_EA=='',])), sep='\n')

##
sq=grepl('sq',qq$hist)
ad=grepl('ad',qq$hist)

##
sq_AA=grepl('sq',qq$hist) & grepl('AA',qq$race)
sq_EA=grepl('sq',qq$hist) & grepl('EA',qq$race)
ad_AA=grepl('ad',qq$hist) & grepl('AA',qq$race)
ad_EA=grepl('ad',qq$hist) & grepl('EA',qq$race)

sq_AA_chtp=data.frame(sort(table(unlist(sapply(qq$chtp[sq_AA], function(x) unlist(strsplit(gsub('[\\(\\) ]','',as.character(x)), '\\,'))  ))), decreasing=T))
colnames(sq_AA_chtp)[1] = 'Chr'
sq_EA_chtp=data.frame(sort(table(unlist(sapply(qq$chtp[sq_EA], function(x) unlist(strsplit(gsub('[\\(\\) ]','',as.character(x)), '\\,'))  ))), decreasing=T))
colnames(sq_EA_chtp)[1] = 'Chr'
ad_AA_chtp=data.frame(sort(table(unlist(sapply(qq$chtp[ad_AA], function(x) unlist(strsplit(gsub('[\\(\\) ]','',as.character(x)), '\\,'))  ))), decreasing=T))
colnames(ad_AA_chtp)[1] = 'Chr'
ad_EA_chtp=data.frame(sort(table(unlist(sapply(qq$chtp[ad_EA], function(x) unlist(strsplit(gsub('[\\(\\) ]','',as.character(x)), '\\,'))  ))), decreasing=T))
colnames(ad_EA_chtp)[1] = 'Chr'


sq_AA_chtp=sq_AA_chtp[as.numeric(as.character(sq_AA_chtp$Var1))<23,]
sq_EA_chtp=sq_EA_chtp[as.numeric(as.character(sq_EA_chtp$Var1))<23,]
ad_AA_chtp=ad_AA_chtp[as.numeric(as.character(ad_AA_chtp$Var1))<23,]
ad_EA_chtp=ad_EA_chtp[as.numeric(as.character(ad_EA_chtp$Var1))<23,]


##
p1<-ggplot(sq_AA_chtp, aes(x=Chr, y=Freq))+geom_bar(stat='identity')+ggtitle('Squamous cell carcinoma of AA')+ylim(x=c(0,max(ad_EA_chtp$Freq, ad_AA_chtp$Freq, sq_EA_chtp$Freq, sq_AA_chtp$Freq) ))
q1<-ggplot(ad_AA_chtp, aes(x=Chr, y=Freq))+geom_bar(stat='identity')+ggtitle('AA-Adenocarcinoma of AA')+ylim(x=c(0,max(ad_EA_chtp$Freq, ad_AA_chtp$Freq, sq_EA_chtp$Freq, sq_AA_chtp$Freq) ))

p2<-ggplot(sq_EA_chtp, aes(x=Chr, y=Freq))+geom_bar(stat='identity')+ggtitle('Squamous cell carcinoma of EA')+ylim(x=c(0,max(ad_EA_chtp$Freq, ad_AA_chtp$Freq, sq_EA_chtp$Freq, sq_AA_chtp$Freq) ))

q2<-ggplot(ad_EA_chtp, aes(x=Chr, y=Freq))+geom_bar(stat='identity')+ggtitle('Adenocarcinoma of EA')+ylim(x=c(0,max(ad_EA_chtp$Freq, ad_AA_chtp$Freq, sq_EA_chtp$Freq, sq_AA_chtp$Freq) ))

##
tiff('/home/sinhas8/Projects/Project_Chromotrypsis/Results_New/chtp_Dist.tiff')
grid.arrange(p1,q1, p2, q2)
dev.off()


##
#qq$chtp=qq$chtp[qq$chtp!==24]
fit1<-survfit(Surv(survival, lungcancer_death_all_years) ~ chtp_Quan, qq[sq_AA,])
fit2<-survfit(Surv(survival, lungcancer_death_all_years) ~ chtp_Quan, qq[sq_EA,])
fit3<-survfit(Surv(survival, lungcancer_death_all_years) ~ chtp_Quan, qq[ad_EA,])
fit4<-survfit(Surv(survival, lungcancer_death_all_years) ~ chtp_Quan, qq[ad_EA,])

p1<-ggsurvplot(fit1, data = qq[sq_AA,], pval=TRUE, title='SCC of AA')
p2<-ggsurvplot(fit2, data = qq[sq_EA,], pval=TRUE, title='SCC of EA')
p3<-ggsurvplot(fit3, data = qq[ad_AA,], pval=TRUE, title='ADC of AA')
p4<-ggsurvplot(fit4, data = qq[ad_EA,], pval=TRUE, title='ADC of EA')

tiff('/home/sinhas8/Projects/Project_Chromotrypsis/Results_New/Kaplan_Mier_CHTP.tiff')
arrange_ggsurvplots(list(p1, p3, p2, p4), ncol=2, nrow=2)
dev.off()

