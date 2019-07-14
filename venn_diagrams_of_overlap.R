########################
#Figure 2B
########################
require(VennDiagram); install.packages('metafolio')
B1=draw.pairwise.venn(area1= LUAD_ovp[[2]][1],
                      area2= LUAD_ovp[[2]][2],
                      cross.area= LUAD_ovp[[2]][3],
                      #                  category= c("AA", "EA"),
                      fill= c('#e41a1c','#377eb8'),
                      main='Amp',
                      cex= 2,
                      scaled=FALSE,
                      alpha=0.8,
                      cat.cex = 2,
                      rotation.degree=180
)
B2= draw.pairwise.venn(area1= LUAD_ovp[[3]][1],
                       area2= LUAD_ovp[[3]][2],
                       cross.area= LUAD_ovp[[3]][3],
                       #                    category= c("AA", "EA"),
                       fill= c('#e41a1c','#377eb8'),
                       main='Recurrent',
                       cex= 2,
                       scaled=FALSE,
                       alpha=0.8,
                       cat.cex = 2, 
                       rotation.degree=180
)
B3= draw.pairwise.venn(area1= LUSC_ovp[[2]][1],
                       area2= LUSC_ovp[[2]][2],
                       cross.area= LUSC_ovp[[2]][3],
                       #category= c("AA", "EA"),
                       fill= c('#e41a1c','#377eb8'),
                       main='Recurrent',
                       cex= 2,
                       scaled=FALSE,
                       alpha=0.8,
                       cat.cex = 2,
                       rotation.degree=180
)
B4= draw.pairwise.venn(area1= LUSC_ovp[[3]][1],
                       area2= LUSC_ovp[[3]][2],
                       cross.area= LUSC_ovp[[3]][3],
                       #category= c("AA", "EA"),
                       fill= c('#e41a1c','#377eb8'),
                       main='Recurrent',
                       cex= 2,
                       scaled=FALSE,
                       alpha=0.8,
                       cat.cex = 2,
                       rotation.degree=180
)
require(cowplot); require(gridExtra)
Figure3_LUAD = plot_grid(grid.arrange(gTree(children=B1), top=textGrob("Amp", gp=gpar(fontsize=20,font=8))),
                         grid.arrange(gTree(children=B2), top=textGrob("Del", gp=gpar(fontsize=20,font=8))),
                         nrow=1)
Figure3_LUSC = plot_grid(grid.arrange(gTree(children=B3), top=textGrob("Amp", gp=gpar(fontsize=20,font=8))),
                         grid.arrange(gTree(children=B4), top=textGrob("Del", gp=gpar(fontsize=20,font=8))),
                         nrow=1)
getwd()
tiff('prep_final_figures/Fig3LUAD_21May.tif', width=300, height = 150)
plot(Figure3_LUAD)
dev.off()
tiff('prep_final_figures/Fig3LUSC_21May.tif', width=300, height = 150)
plot(Figure3_LUSC)
dev.off()
