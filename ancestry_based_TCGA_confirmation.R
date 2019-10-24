################################################################################
####Process the inferred ancestry
################################################################################
setwd('/Users/sinhas8/Project_Chromotrypsis/2.Data/inferred_ancestry/')
IA=lapply(1:length(list.files()), function(x) readRDS(list.files()[x]))
IA=do.call(rbind, IA)
tcga=read.csv('/Users/sinhas8/Project_Chromotrypsis/Results_New/Nov_28/Supp_Table3.csv')
IA_withgi=cbind(IA, tcga[match(IA$sample, tcga$info_tcga.Sample_Name),])
tcga=IA_withgi
table(as.character(tcga$info_tcga.hist))
tcga=tcga[!is.na(tcga$info_tcga.hist),]

tcga_list=split(tcga, as.character(tcga$info_tcga.hist))
head(tcga_list[[1]]$inferred_ancestry== as.character(tcga_list[[1]]$race))
correct_ancestry=sapply(tcga_list, function(x) !sum(as.character(tcga_list[[1]]$race) == x$inferred_ancestry)>(0.5*nrow(x)) )
corr_ancestry=do.call(rbind, lapply(tcga_list[correct_ancestry], function(x) data.frame(x$file_id, factor(x$inferred_ancestry, labels = c('EA','AA'))) ))
table(tcga$inferred_ancestry[match(corr_ancestry$x.file_id, tcga$file_id)], corr_ancestry$factor.x.inferred_ancestry..labels...c..EA....AA...)
tcga$inferred_ancestry[match(corr_ancestry$x.file_id, tcga$file_id)] = corr_ancestry$factor.x.inferred_ancestry..labels...c..EA....AA...


range01 <- function(x){(x-min(x))/(max(x)-min(x))}
scaling_cancerType<-function(quan1=gi, quan2=hist){
  unlist(lapply(split(quan1, quan2), function(x) range01(x)))
}
tcga$Normalized_gi=scaling_cancerType(tcga$CNV_Burden, tcga$info_tcga.hist)
tcga$Normalized_hrd.loh=scaling_cancerType(tcga$HRD_by_LOH, tcga$info_tcga.hist)
tcga$Normalized_hrd.LST=scaling_cancerType(tcga$HRD_by_LST, tcga$info_tcga.hist)
tcga$Normalized_hrd.AIL=scaling_cancerType(tcga$HRD_by_AIL, tcga$info_tcga.hist)
tcga$Normalized_Chromothripsis_Presence=scaling_cancerType(tcga$CHTP_Canonical_Definition_Presence, tcga$info_tcga.hist)
# colnames(tcga)[6]='race'
# levels(tcga$race)[2]=c('EA')
# colnames(tcga)[7]='hist'
##################
###Getting Legend
##################
ylim1=c(0,1)
leg=get_legend(ggplot(tcga, aes(x=race,y=Normalized_gi, fill=inferred_ancestry))+
                 geom_violin(trim=FALSE)+
                 geom_boxplot(width=0.3, fill="white")+
                 labs(title="GI in PanCan (TCGA)",x="inferred_ancestry", y = "Scaled GI",
                      fill='inferred ancestry')+
                 theme_classic(base_size = 25)+
                 stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                                    label.y=ylim1[2], size = 10)+
                 theme(legend.box = "horizontal", legend.position = "bottom",
                       legend.text=element_text(size=20), legend.title=element_text(size=25),
                       legend.key.size = unit(2,"line"))+
                 guides(fill=guide_legend(title="Inferred Ancestry"))+
                 scale_fill_manual(values = c('#e41a1c','#377eb8'), 
                                   labels=c('African Americans', 'European Americans'),
                                   name="Inferred Ancestry")
)

##################
#S1TCGA:GI
##################
Plot1=ggplot(tcga, aes(x=as.character(info_tcga.hist),y=Normalized_gi, fill=inferred_ancestry))+
  geom_boxplot(data=tcga, aes(fill=inferred_ancestry))+
  labs(title="GI in 23 cancer types (TCGA)",x="inferred_ancestry", y = "Scaled GI", 
       fill='inferred ancestry')+
  facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                     #label.y=ylim1[2],
                     size = 10)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)#+
#coord_cartesian(ylim = ylim1*1.15)
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/FigureSupp10.tiff', width = 1800, height = 650)
p1=plot_grid(leg, Plot1, ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()

setwd('/Users/sinhas8/Project_Chromotrypsis/')
tiff('prep_final_figures/FigureSupp10_HRDLST.tif', width = 1800, height = 650)
p2=plot_grid(leg, 
          ggplot(tcga, aes(x=as.character(info_tcga.hist),y=Normalized_hrd.LST, fill=inferred_ancestry))+
            geom_boxplot(data=tcga, aes(fill=inferred_ancestry))+
            labs(title="HRD by LST in 23 cancer types (TCGA)",x="inferred_ancestry", y = "Scaled HRD by LST")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
          , ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()

tiff('prep_final_figures/FigureSupp10_HRDAIL.tif', width = 1800, height = 650)
p3=plot_grid(leg, 
          ggplot(tcga, aes(x=as.character(info_tcga.hist),y=Normalized_hrd.AIL, fill=inferred_ancestry))+
            geom_boxplot(data=tcga, aes(fill=inferred_ancestry))+
            labs(title="HRD by AIL in 23 cancer types (TCGA)",x="inferred_ancestry", y = "Scaled HRD by AIL")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
          , ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()

tiff('prep_final_figures/FigureSupp10_HRDLOH.tif', width = 1800, height = 650)
p4=plot_grid(leg, 
          ggplot(tcga, aes(x=as.character(info_tcga.hist),y=Normalized_hrd.loh, fill=inferred_ancestry))+
            geom_boxplot(data=tcga, aes(fill=inferred_ancestry))+
            labs(title="HRD by LOH in 23 cancer types (TCGA)",x="inferred_ancestry", y = "Scaled HRD by LOH")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
          , ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
tiff('prep_final_figures/FigureSupp10_HRDnetsum.tif', width = 1800, height = 650)
p5=plot_grid(leg, 
          ggplot(tcga, aes(x=as.character(info_tcga.hist),
                           y=range01(Normalized_hrd.loh+Normalized_hrd.AIL+Normalized_hrd.LST),
                           fill=inferred_ancestry))+
            geom_boxplot(data=tcga, aes(fill=inferred_ancestry))+
            labs(title="HRD by netsum of 3 measure in 23 cancer types (TCGA)",x="inferred_ancestry", y = "Scaled HRD by all three sum")+
            facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
            theme_classic(base_size = 25)+
            stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                               label.y=ylim1[2], size = 10)+
            guides(fill=FALSE)+
            scale_fill_brewer(palette=colorType_Set)
          , ncol = 1, rel_heights = c(1/20, 19/20))
dev.off()

chtp_mat=aggregate(CHTP_Canonical_Definition_Presence ~ inferred_ancestry+info_tcga.hist+info_tcga.CellofOrigin+info_tcga.Tissue_Type,
                   function(x) sum(x)/length(x), 
                   data = tcga)
chtp_mat=chtp_mat[unlist(lapply(seq(1, nrow(chtp_mat), by = 2), function(x)
  rep(sum(chtp_mat$CHTP_Canonical_Definition_Presence[x:x+1])>0, 2) )),]
tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/FigureSupp10_CHTP.tif', 
     width = 1350, height = 450)
p6=plot_grid(leg,
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

tiff('/Users/sinhas8/Project_Chromotrypsis/prep_final_figures/FigureSupp10.tif', width = 2500, height = 1500)
plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3, ncol=2, labels = 'AUTO', label_size = 30)
dev.off()


