##################
#S1TCGA:Plotting Supplemntary
##################
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
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                     #label.y=ylim1[2],
                     size = 10)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)#+
  #coord_cartesian(ylim = ylim1*1.15)
tiff('prep_final_figures/FigureSupp1.tif', width = 1800, height = 650)
plot_grid(leg, Plot1, ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()
##################
#S2TCGA:GI gain and loss
##################
df=read.csv('/Users/sinhas8/df_GIlossandgain.csv')
df$Normalized_gi=scaling_cancerType(df$gi, df$hist)
df$Normalized_gi_gain=scaling_cancerType(df$GI_gain, df$hist)

Plot2=ggplot(df, aes(x=as.character(hist),y=Normalized_gi, fill=race))+
  geom_boxplot(data=df, aes(fill=race))+
  labs(title="loss- GI in 23 cancer types (TCGA)",x="Race", y = "Scaled loss- GI")+
  facet_grid( ~Tissue_Type + CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                     label.y=ylim1[2], size = 10)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)#+

tiff('prep_final_figures/FigureSupp2A.tif', width = 1800, height = 650)
plot_grid(leg, Plot2, ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()

Plot3=ggplot(df, aes(x=as.character(hist),y=Normalized_gi_gain, fill=race))+
  geom_boxplot(data=df, aes(fill=race))+
  labs(title="gain- GI in 23 cancer types (TCGA)",x="Race", y = "Scaled gain- GI")+
  facet_grid( ~Tissue_Type + CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                     label.y=ylim1[2], size = 10)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)#+

tiff('prep_final_figures/FigureSupp2B.tif', width = 1800, height = 650)
plot_grid(leg, Plot3, ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()
##################
#S3TCGA:HRD across cancer types
##################
# df$Normalized_Chromothripsis_Presence=scaling_cancerType(df$CHTP_Canonical_Definition_Presence, df$hist)
tiff('prep_final_figures/FigureSupp3A.tif', width = 1800, height = 650)
plot_grid(leg, 
          ggplot(tcga, aes(x=as.character(hist),y=Normalized_hrd.LST, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="HRD by LST in 23 cancer types (TCGA)",x="Race", y = "Scaled HRD by LST")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
          , ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()
tiff('prep_final_figures/FigureSupp3B.tif', width = 1800, height = 650)
plot_grid(leg, 
          ggplot(tcga, aes(x=as.character(hist),y=Normalized_hrd.AIL, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="HRD by AIL in 23 cancer types (TCGA)",x="Race", y = "Scaled HRD by AIL")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
          , ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()

tiff('prep_final_figures/FigureSupp3C.tif', width = 1800, height = 650)
plot_grid(leg, 
          ggplot(tcga, aes(x=as.character(hist),y=Normalized_hrd.loh, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="HRD by LOH in 23 cancer types (TCGA)",x="Race", y = "Scaled HRD by LOH")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
          , ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
tiff('prep_final_figures/FigureSupp3D.tif', width = 1800, height = 650)
plot_grid(leg, 
          ggplot(tcga, aes(x=as.character(hist),
                           y=range01(Normalized_hrd.loh+Normalized_hrd.AIL+Normalized_hrd.LST),
                           fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="HRD by netsum of 3 measure in 23 cancer types (TCGA)",x="Race", y = "Scaled HRD by all three sum")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
          , ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()
##################
#S4TCGA:HRD in PanCancer
##################
# df$Normalized_Chromothripsis_Presence=scaling_cancerType(df$CHTP_Canonical_Definition_Presence, df$hist)
tiff('prep_final_figures/FigureSupp4A.tif', width = 400, height = 400)
ggplot(tcga, aes(x=race,y=Normalized_hrd.LST, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="",x="Race", y = "Scaled HRD by LST")+
            #facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.2,
                               method.args = list(alternative = "greater"),
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
dev.off()
tiff('prep_final_figures/FigureSupp4B.tif', width = 400, height = 400)
ggplot(tcga, aes(x=race,
                           y=Normalized_hrd.AIL, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="",x="Race", y = "Scaled HRD by AIL")+
            #facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.2,
                     method.args = list(alternative = "greater"),
                     label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
dev.off()

tiff('prep_final_figures/FigureSupp4C.tif', width = 400, height = 400)
          ggplot(tcga, aes(x=race, y=Normalized_hrd.loh, fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="",x="Race", y = "Scaled HRD by LOH")+
            #facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.2,
                               method.args = list(alternative = "greater"),
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
dev.off()

tiff('prep_final_figures/FigureSupp4D.tif', width = 400, height = 400)
          ggplot(tcga, aes(x=race,
                           y=range01(Normalized_hrd.loh+Normalized_hrd.AIL+Normalized_hrd.LST),
                           fill=race))+
            geom_boxplot(data=tcga, aes(fill=race))+
            labs(title="",x="Race", y = "Scaled HRD by all three sum")+
            #facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.2,
                               method.args = list(alternative = "greater"),
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
dev.off()


##################
#HRD by PanCacer
##################
TCGAsig3=TCGAsig[TCGAsig$Signature=='Signature.3',]
TCGAsig3=TCGAsig3[grep('TCGA',TCGAsig3$Project_Name),]
TCGAsig3$hist=sapply(TCGAsig3$project_code, 
                     function(x) strsplit(as.character(x), split = '-')[[1]][1])
df=cbind(TCGAsig3, tcga[match(TCGAsig3$hist,tcga$hist),
                        c('info_tcga.Tissue_Type','info_tcga.CellofOrigin')])
head(df)
df=df[!is.na(df$info_tcga.Tissue_Type),]
tiff('prep_final_figures/FigureSupp3E.tif', width = 2000, height = 650)
plot_grid(leg, 
          ggplot(df, aes(x=as.character(hist),y=Contribution, fill=race))+
            geom_boxplot(data=df, aes(fill=race))+
            labs(title="HRD by Signature",x="Race", y = "Scaled HRD by Sig3")+
            theme_classic(base_size = 25)+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                               method.args = list(alternative = "greater"),
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
          , ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()


##################
#S4EPanCancer Sig3 Quan
##################
tiff('prep_final_figures/FigureSupp4E.tif', width = 400, height = 400)
ggplot(df, aes(x=race,y=Contribution, fill=race))+
  geom_boxplot(data=df, aes(fill=race))+
  labs(title="HRD by Signature",x="Race", y = "Scaled HRD by Sig3")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                     method.args = list(alternative = "greater"),
                     label.y=0.15, size = 10)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)+
  ylim(c(0, 0.2))
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
#NCIMD- CHTP in lung cancer
##################
mat$chtp_Quan=as.numeric(mat$chtp_Quan>0)
fisher.test(cbind(table(mat$chtp_Quan[mat$hist=='LUAD' & mat$race=='EA']),
                  table(mat$chtp_Quan[mat$hist=='LUAD' & mat$race=='AA'])),
            alternative = 'g')
chtp_mat=aggregate(chtp_Quan ~ race, function(x) sum(x)/length(x), data = mat[mat$hist=='LUAD',])
F2=ggplotGrob(ggplot(chtp_mat, aes(y=chtp_Quan, x=race, fill=race))+
                geom_bar(stat = 'identity')+
                labs(title="CHTP in LUAD (NCI-MD) ",x="Race", y = "Proportion with Chromothripsis")+
                theme_classic(base_size = 25)+
                guides(fill=FALSE)+
                annotate("text", x = 1.5, y = max(chtp_mat[,2])*1.05, label = "p = 0.12", size=10)+
                coord_cartesian(ylim = c(0, max(chtp_mat[,2])*1.15))+
                scale_fill_brewer(palette=colorType_Set) 
)

tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/SuppFigS5.tif', 
     width = 1350, height = 450)
plot_grid(D2, E2, F2, align='h', nrow=1,
                    labels = c(LETTERS[1:3]), label_size = 35 ) 
dev.off()

##################
#TCGA- CHTP across cancer types
##################
fisher.test(cbind(table(tcga$CHTP_Canonical_Definition_Presence[tcga$race=='AA' & lung_can]),
                  table(tcga$CHTP_Canonical_Definition_Presence[tcga$race=='EA' & lung_can])),
            alternative = 'l')
chtp_mat=aggregate(CHTP_Canonical_Definition_Presence ~ inferred_ancestry+info_tcga.hist+info_tcga.CellofOrigin+info_tcga.Tissue_Type,
                   function(x) sum(x)/length(x), 
                   data = tcga)
chtp_mat=chtp_mat[unlist(lapply(seq(1, nrow(chtp_mat), by = 2), function(x)
  rep(sum(chtp_mat$CHTP_Canonical_Definition_Presence[x:x+1])>0, 2) )),]
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/SuppFigS6.tif', 
     width = 1350, height = 450)
plot_grid(leg,
          ggplotGrob(ggplot(chtp_mat, aes(y=CHTP_Canonical_Definition_Presence, x=info_tcga.hist, fill=inferred_ancestry))+
                geom_bar(stat = 'identity', position="dodge")+
                labs(title="CHTP across cancer types from TCGA",x="Race", y = "Proportion with Chromothripsis")+
                theme_classic(base_size = 25)+
                facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
                guides(fill=FALSE)+
                #annotate("text", x = 1.5, y = max(chtp_mat[,2])*1.05, label = "p = 0.11", size=10)+
                #coord_cartesian(ylim = c(0, max(chtp_mat[,2])*1.15))+
                scale_fill_brewer(palette=colorType_Set)),
          ncol = 1, rel_heights = c(1/20, 19/20)
          )          
dev.off()

##################
#TCGA- CHTP across chromosomes
##################
