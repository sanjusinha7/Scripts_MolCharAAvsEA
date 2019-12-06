# Figure 3
########################
# Venn Plots
########################
require(VennDiagram)
B1=draw.pairwise.venn(area1= 30,
                      area2= 27,
                      cross.area= 12,
                      #                  category= c("AA", "EA"),
                      fill= c('#e41a1c','#377eb8'),
                      main='Recurrent',
                      cex= 2,
                      scaled=FALSE,
                      alpha=0.8,
                      cat.cex = 2
)
B2= draw.pairwise.venn(area1= 45,
                       area2= 54,
                       cross.area= 24,
                       #                    category= c("AA", "EA"),
                       fill= c('#e41a1c','#377eb8'),
                       main='Recurrent',
                       cex= 2,
                       scaled=FALSE,
                       alpha=0.8,
                       cat.cex = 2, 
                       rotation.degree=180
)
B3= draw.pairwise.venn(area1= 18,
                       area2= 34,
                       cross.area= 10,
                       #category= c("AA", "EA"),
                       fill= c('#e41a1c','#377eb8'),
                       main='Recurrent',
                       cex= 2,
                       scaled=FALSE,
                       alpha=0.8,
                       cat.cex = 2,
                       rotation.degree=180
)
B4= draw.pairwise.venn(area1= 32,
                       area2= 34,
                       cross.area= 12,
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
Figure3_LUAD = 
  Figure3_LUSC = plot_grid(grid.arrange(gTree(children=B3), top=textGrob("Amp", gp=gpar(fontsize=20,font=8))),
                           grid.arrange(gTree(children=B4), top=textGrob("Del", gp=gpar(fontsize=20,font=8))),
                           nrow=1)
getwd()

pdf('../prep_final_figures/Fig3LUAD.pdf', width=4, height = 2)
plot(Figure3_LUAD)
dev.off()

pdf('../prep_final_figures/Fig3LUSC.pdf', width=4, height = 2)
plot(Figure3_LUSC)
dev.off()


names(postscriptFonts())