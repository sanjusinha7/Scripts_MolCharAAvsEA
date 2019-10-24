#######Write Supplementary table 14(initial) - table with population distributions
df2write=data.frame(Cancer_Type=as.character(aggregate(X ~ hist, data=tcga, function(x) length(x) )[,1]),
           Total_Samples=aggregate(X ~ hist, data=tcga, function(x) length(x) )[,2],
           AA_tumors=aggregate(X ~ hist+race, data=tcga, function(x) length(x) )[1:23,3],
           EA_tumors=aggregate(X ~ hist+race, data=tcga, function(x) length(x) )[24:46,3])

write.csv(df2write, './2.Data/Supptable14.csv')
df2write
###Add Median GI, HRD and CHTP frequency.
custom_agg<-function(agg1){
  do.call(cbind, split(agg1[,3], agg1[,2]))
}

df2write=cbind(df2write, data.frame(Median_norm.GI=custom_agg(aggregate(tcga$Normalized_gi ~ hist+race, data=tcga, function(x) median(x) )),
                                    Median_norm.hrd_loh=custom_agg(aggregate(Normalized_hrd.loh ~ hist+race, data=tcga, function(x) median(x) )),
                                    Median_norm.hrd_lst=custom_agg(aggregate(tcga$Normalized_hrd.LST ~ hist+race, data=tcga, function(x) median(x) )),
                                    Median_norm.hrd_ail=custom_agg(aggregate(tcga$Normalized_hrd.AIL ~ hist+race, data=tcga, function(x) median(x) )) ))

head(df2write)