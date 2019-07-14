###Pre-Processing for HRD-pathway

Mut_raw=read.csv('/home/sinhas/Downloads/portal-mutation-2018-05-09.csv', sep=',')
Mut_Filtered = Mut_raw[which(!(Mut_raw$variant_classification=='Silent')),]
Mut_FF=Mut_Filtered[,c(1,15)]
qq=do.call(rbind,lapply(split(Mut_FF[,1], Mut_FF[,2]), table))
mat=t(qq)


