##For OncoScan 
Onco.Peaks_Amp_EA = read.csv('/home/sinhas8/Both_EA/amp_genes.conf_90.txt', fill=T, sep='\t', header=F)
Onco.Peaks_Del_EA = read.csv('/home/sinhas8/Both_EA/del_genes.conf_90.txt', fill=T, sep='\t', header=F)

Onco.Peaks_Amp_AA = read.csv('/home/sinhas8/Both_AA/amp_genes.conf_90.txt', fill=T, sep='\t', header=F)
Onco.Peaks_Del_AA = read.csv('/home/sinhas8/Both_AA/del_genes.conf_90.txt', fill=T, sep='\t', header=F)

#Genes
Gene_q.Score<-function(Genes_OncoScan){
	q_values=Genes_OncoScan[2,-1]
	tt=Genes_OncoScan[-c(1:4),-1]
	GeneScore = mapply(function(a, b) data.frame(GeneName=unique(a[a!='']),Q_value=b) ,tt, q_values, SIMPLIFY=F)
	GS_df=do.call(rbind, GeneScore)
}

#Extracting the genes::
Onco_Amp_AA=Gene_q.Score(Onco.Peaks_Amp_AA)
Onco_Del_AA=Gene_q.Score(Onco.Peaks_Del_AA)
write.csv(Onco_Amp_AA$GeneName, '/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/Onco_Amp_AA.csv')
write.csv(Onco_Del_AA$GeneName, '/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/Onco_Del_AA.csv')

Onco_Amp_EA=Gene_q.Score(Onco.Peaks_Amp_EA)
Onco_Del_EA=Gene_q.Score(Onco.Peaks_Del_EA)
write.csv(Onco_Amp_EA$GeneName, '/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/Onco_Amp_EA.csv')
write.csv(Onco_Del_EA$GeneName, '/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/Onco_Del_EA.csv')

Amp_inboth=data.frame(Onco_Amp_EA$GeneName[na.omit(match(Onco_Amp_AA$GeneName, Onco_Amp_EA$GeneName))])
Del_inboth=data.frame(Onco_Del_EA$GeneName[na.omit(match(Onco_Del_AA$GeneName, Onco_Del_EA$GeneName))])
write.csv(Amp_inboth, '/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/Amp_inboth.csv')
write.csv(Del_inboth, '/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/Del_inboth.csv')


Onco_Amp_total=unique(unlist(list(unlist(Onco_Amp_AA$GeneName), unlist(Onco_Amp_EA$GeneName))))
Onco_Del_total=unique(unlist(list(unlist(Onco_Del_AA$GeneName), unlist(Onco_Del_EA$GeneName))))
write.csv(Onco_Amp_total, '/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/Onco_Amp_total.csv')
write.csv(Onco_Del_total, '/home/sinhas8/Projects/Project_Chromotrypsis/2.Data/Onco_Del_total.csv')



