##################
#S1TCGA:Plotting Supplemntary
##################
###Add Sample Size
levels(df$hist)=
  apply(aggregate(X ~ hist, data = df, length), 1, function(x) paste0(x[1], '\n','(',x[2],')',sep=''))
tcga$hist_wdoutSampleSize=tcga$hist
levels(tcga$hist)=
  apply(aggregate(X ~ hist, data = tcga, length), 1, function(x) paste0(x[1], '\n','(',x[2],')',sep=''))

##################
###Getting Legend
##################
ylim1=c(0,1)
leg=get_legend(ggplot(tcga, aes(x=race,y=Normalized_gi, fill=race))+
                 geom_violin(trim=FALSE)+
                 geom_boxplot(width=0.3, fill="white")+
                 labs(title="GI in PanCan (TCGA)",x="Race", y = "Scaled GI")+
                 theme_classic(base_size = 25)+
                 stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                                    label.y=ylim1[2], size = 10)+
                 theme(legend.box = "horizontal", legend.position = "bottom",
                       legend.text=element_text(size=20), legend.title=element_text(size=25),
                       legend.key.size = unit(2,"line"))+
                 guides(fill=guide_legend(title="Race"))+
                 scale_fill_manual(values = c('#e41a1c','#377eb8'), 
                                   labels=c('African Americans', 'European Americans'),
                                   name="Race")
)

##################
#S1TCGA:GI
##################
Plot1=ggplot(tcga, aes(x=as.character(hist),y=Normalized_gi, fill=race))+
  geom_boxplot(data=tcga, aes(fill=race))+
  labs(title="GI in 23 cancer types (TCGA)",x="Race", y = "Scaled GI")+
  facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     label.y=c(1.03, 1.1, 1.03, 1.1, 1.03, 1.1,1.03),
                     size = 6)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)#+

tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Extended_Figure5_Nov5.tif',
     width = 1800, height = 650)
plot_grid(leg, Plot1, ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()
##################
#S2TCGA:GI gain and loss
##################
df=read.csv('/Users/sinhas8/df_GIlossandgain.csv')
df$Normalized_gi=scaling_cancerType(df$gi, df$hist)
df$Normalized_gi_gain=scaling_cancerType(df$GI_gain, df$hist)

Plot2=ggplot(df_wdSampleSize, aes(x=as.character(hist),y=Normalized_gi, fill=race))+
  geom_boxplot()+
  labs(title="loss- GI in 23 cancer types (TCGA)",x="Race", y = "Scaled loss- GI")+
  facet_grid( ~Tissue_Type + CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     label.y=c(1.03, 1.1, 1.03, 1.1, 1.03, 1.1,1.03),
                     size = 5.5)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)#+

Plot3=ggplot(df_wdSampleSize, aes(x=hist,y=Normalized_gi_gain, fill=race))+
  geom_boxplot()+
  labs(title="gain- GI in 23 cancer types (TCGA)",x="Race", y = "Scaled gain- GI")+
  facet_grid( ~Tissue_Type + CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     label.y=c(1.03, 1.1, 1.03, 1.1, 1.03, 1.1,1.03),
                     size = 5.5)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)

tiff('prep_final_figures/Ext_Fig6.tif', width = 1800, height = 1300)
plot_grid(leg, Plot2, Plot3, ncol = 1, rel_heights = c(1/20, 9.5/20, 9.5/20), 
          labels = c('','A', 'B'), label_size = 30)
dev.off()
##################
#S3
##################
Fig3A=ggplot(tcga, aes(x=as.character(hist),y=Normalized_hrd.LST, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="HRD by LST in 23 cancer types (TCGA)",x="Race", y = "Scaled HRD by LST")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                               label.y=c(1.03, 1.1, 1.03, 1.1, 1.03, 1.1,1.03),
                               size = 6)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
# dev.off()
# tiff('prep_final_figures/FigureSupp3B.tif', width = 1800, height = 650)
Fig3B=ggplot(tcga, aes(x=as.character(hist),y=Normalized_hrd.AIL, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="HRD by AIL in 23 cancer types (TCGA)",x="Race", y = "Scaled HRD by AIL")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                               label.y=c(1.03, 1.1, 1.03, 1.1, 1.03, 1.1,1.03),
                               size = 6)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
# dev.off()

# tiff('prep_final_figures/FigureSupp3C.tif', width = 1800, height = 650)
Fig3C=ggplot(tcga, aes(x=as.character(hist),y=Normalized_hrd.loh, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="HRD by LOH in 23 cancer types (TCGA)",x="Race", y = "Scaled HRD by LOH")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                               label.y=c(1.03, 1.1, 1.03, 1.1, 1.03, 1.1,1.03),
                               size = 6)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
# dev.off()
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# tiff('prep_final_figures/FigureSupp3D.tif', width = 1800, height = 650)
Fig3D=ggplot(tcga, aes(x=as.character(hist),
                           y=range01(Normalized_hrd.loh+Normalized_hrd.AIL+Normalized_hrd.LST),
                           fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="HRD by netsum of 3 measure in 23 cancer types (TCGA)",x="Race", y = "Scaled HRD by all three sum")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                               label.y=c(1.03, 1.1, 1.03, 1.1, 1.03, 1.1,1.03),
                               size = 6)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)

#Load Signature information
TCGAsig=read.csv('/Users/sinhas8/Downloads/signature_profile_sample.txt', sep='\t')
levels(TCGAsig$race)[c(3,6)]=c('AA','EA')
TCGAsig3=TCGAsig[TCGAsig$Signature=='Signature.3',]
TCGAsig3=TCGAsig3[grep('TCGA',TCGAsig3$Project_Name),]
TCGAsig3$hist=sapply(TCGAsig3$project_code, 
                     function(x) strsplit(as.character(x), split = '-')[[1]][1])
df=cbind(TCGAsig3, tcga[match(TCGAsig3$hist,tcga$hist_wdoutSampleSize),
                        c('info_tcga.Tissue_Type','info_tcga.CellofOrigin', 'hist_wdSampleSize')])
df=df[!is.na(df$info_tcga.Tissue_Type),]

Fig3E=ggplot(df, aes(x=as.character(hist_wdSampleSize),y=Contribution, fill=race))+
  geom_boxplot(data=df, aes(fill=race))+
  labs(title="HRD by Signature",x="Race", y = "Scaled HRD by Sig3")+
  theme_classic(base_size = 25)+
  facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin,
              scales = "free", space = "free_x")+
  stat_compare_means(method = "wilcox.test",
                     label = "p", 
                     label.x = 1.5,
                     label.y=c(1.03, 1.1, 1.03, 1.1, 1.03, 1.1,1.03),
                     method.args = list(alternative = "greater"),
                     size = 6)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)


tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Extended_Figure8_Nov5.tif',
     width = 3400, height = 1900)
plot_grid(leg, plot_grid(Fig3A, Fig3B, Fig3C,
                         Fig3D, Fig3E, nrow = 3,
                         ncol=2, labels = 'AUTO', label_size = 30),
          nrow=2,
          rel_heights = c(1/50, 49/50))
dev.off()

##################
#S4
##################
Fig4A=ggplot(tcga, aes(x=race,y=Normalized_hrd.LST, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="",x="Race", y = "Scaled HRD \nby LST")+
            #facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                               label.y=c(1.03, 1.1, 1.03, 1.1, 1.03, 1.1,1.03),
                     method.args = list(alternative = "greater"),
                     size = 6)+
              guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
# dev.off()
# tiff('prep_final_figures/FigureSupp4B.tif', width = 400, height = 400)
Fig4B=ggplot(tcga, aes(x=race,
                           y=Normalized_hrd.AIL, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="",x="Race", y = "Scaled HRD\n by AIL")+
            #facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                               label.y=c(1.03, 1.1, 1.03, 1.1, 1.03, 1.1,1.03),
                     method.args = list(alternative = "greater"),
                     size = 6)+
              guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
# dev.off()

# tiff('prep_final_figures/FigureSupp4C.tif', width = 400, height = 400)
Fig4C=          ggplot(tcga, aes(x=race, y=Normalized_hrd.loh, fill=race))+
  geom_boxplot(data=tcga, aes(fill=race))+
  labs(title="",x="Race", y = "Scaled HRD \nby LOH")+
  #facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.2,
                     method.args = list(alternative = "greater"),
                     label.y=ylim1[2],
                     size = 6)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)
# dev.off()

# tiff('prep_final_figures/FigureSupp4D.tif', width = 400, height = 400)
Fig4D=          ggplot(tcga, aes(x=race,
                                 y=range01(Normalized_hrd.loh+Normalized_hrd.AIL+Normalized_hrd.LST),
                                 fill=race))+
  geom_boxplot(data=tcga, aes(fill=race))+
  labs(title="",x="Race", y = "Scaled HRD by\nall three sum")+
  #facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.2,
                     method.args = list(alternative = "greater"),
                     label.y=ylim1[2],
                     size = 6)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)

Fig4E=ggplot(df, aes(x=race,y=Contribution, fill=race))+
  geom_boxplot(data=df, aes(fill=race))+
  labs(x="Race", y = "Scaled HRD \nby Sig3")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.2,
                     method.args = list(alternative = "greater"),
                     size = 6)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)+
  ylim(c(0, 0.2))

table(df$race)
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/Extended_Figure7_Nov5.tif',
     width = 800, height = 600)
plot_grid(leg, plot_grid(Fig4A, Fig4B, Fig4C,
                         Fig4D, Fig4E, nrow = 2,
                         ncol=3, labels = 'AUTO',  label_size = 30),
          nrow=2,
          rel_heights = c(3/50, 47/50))
dev.off()

##################
#LUAD GI, HRD and CHTP
##################
##################
#NCIMDProcessing
##################
mat=read.csv('/users/sinhas8/Project_Chromotrypsis/2.Data/Corrected_NCIMD_HRD_by_LOH_and_GI.csv')
mat=mat[mat$hist=='adeno' | mat$hist =='sq',]
mat$hist=as.character(mat$hist)
mat$hist=factor(mat$hist)
levels(mat$hist)=c('LUAD', 'LUSC')
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
scaling_cancerType<-function(quan1=gi, quan2=hist){
  unlist(lapply(split(quan1, quan2), function(x) range01(x)))
}
mat=mat[order(mat$hist),]
mat$Normalized_gi=scaling_cancerType(mat$CNV_burden, mat$hist)
mat$Normalized_hrd.loh=scaling_cancerType(mat$HRD_by_LOH, mat$hist)
mat$Normalized_Chromothripsis_Presence=scaling_cancerType(mat$chtp_Quan, mat$hist)

##################
#NCIMD- GI in lung cancer
##################
ylim1=c(0,1)
D2=ggplotGrob(ggplot(mat[mat$hist=='LUAD',], aes(x=race, y=Normalized_gi, fill=race))+
                geom_violin(trim=FALSE)+
                geom_boxplot(width=0.3, fill="white")+
                labs(title="GI in LUAD (NCI-MD)",x="Race", y = "Scaled GI")+
                theme_classic(base_size = 25)+
                stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                                   label.y=ylim1[2], size = 10)+
                guides(fill=FALSE)+
                scale_fill_brewer(palette=colorType_Set)+
                coord_cartesian(ylim = ylim1*1.15)
)
##################
#NCIMD- GI in lung cancer
##################
E2=ggplotGrob(ggplot(mat[mat$hist=='LUAD',], aes(x=race, y=Normalized_hrd.loh, fill=race))+
                geom_violin(trim=FALSE)+
                geom_boxplot(width=0.3, fill="white")+
                labs(title="HRD in LUAD (NCI-MD)",x="Race", y = "Scaled HRD")+
                theme_classic(base_size = 25)+
                stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                                   label.y=ylim1[2], size = 10)+
                guides(fill=FALSE)+
                scale_fill_brewer(palette=colorType_Set)+
                coord_cartesian(ylim = ylim1*1.15))
##################
#TCGA- CHTP across cancer types
##################
##identify p-value 
chtp_mat=aggregate(CHTP_Canonical_Definition_Presence ~ race+info_tcga.hist_wdSampleSize+
                     info_tcga.CellofOrigin+info_tcga.Tissue_Type,
                   function(x) sum(x)/length(x), 
                   data = tcga)
chtp_mat=chtp_mat[unlist(lapply(seq(1, nrow(chtp_mat), by = 2), function(x)
  rep(sum(chtp_mat$CHTP_Canonical_Definition_Presence[x:x+1])>0, 2) )),]

Pvalues=data.frame(sig=sapply(split(tcga, tcga$info_tcga.hist), function(x) 
  err_handle(fisher.test(t(cbind(table(x$CHTP_Canonical_Definition_Presence[x$race=='EA']),
                                 table(x$CHTP_Canonical_Definition_Presence[x$race=='AA']))),
                         alternative = 'g')[1]$p.value)))
sig_df=data.frame(info_tcga.hist_wdSampleSize=levels(tcga$info_tcga.hist_wdSampleSize)[!is.na(Pvalues$sig)],
                  sig=round(na.omit(Pvalues$sig), 2),
                  ystar=c(0.65, 0.70),
                  race='AA',
                  #levels(tcga$info_tcga.hist)[!is.na(Pvalues$sig)]
                  info_tcga.Tissue_Type=tcga$info_tcga.Tissue_Type[
                    match(sig_df$info_tcga.hist_wdSampleSize,
                          tcga$info_tcga.hist_wdSampleSize)],
                  info_tcga.CellofOrigin=tcga$info_tcga.CellofOrigin[
                    match(sig_df$info_tcga.hist_wdSampleSize,
                          tcga$info_tcga.hist_wdSampleSize)]
)
sig_df=sig_df[order(sig_df$info_tcga.Tissue_Type, sig_df$info_tcga.CellofOrigin),]
sig_df$xstar=c(1, 1, 1, 2, 1, 1, 2, 1)
Ext_Fig2A=plot_grid(leg, ggplot(chtp_mat, aes(y=CHTP_Canonical_Definition_Presence, 
                                       x=info_tcga.hist_wdSampleSize, fill=race))+
               geom_bar(stat = 'identity', position="dodge")+
               labs(title="CHTP across cancer types from TCGA",x="Race", y = "Proportion with Chromothripsis")+
               theme_classic(base_size = 25)+
               facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
               guides(fill=FALSE)+
               #annotate("text", x = 1.5, y = max(chtp_mat[,2])*1.05, label = "p = 0.11", size=10)+
               #coord_cartesian(ylim = c(0, max(chtp_mat[,2])*1.15))+
               scale_fill_brewer(palette=colorType_Set)+
               geom_text(aes(label=round(as.numeric(CHTP_Canonical_Definition_Presence), 2)),
                         position=position_dodge(width=0.9), vjust=-0.25, size=6)+
               geom_text(data=sig_df, aes(x=xstar, y=ystar, label=paste0('p=', sig)), size=8),
             ncol = 1, rel_heights = c(1/20, 19/20)
)          

##################
#TCGA- CHTP across chromosomes
##################
