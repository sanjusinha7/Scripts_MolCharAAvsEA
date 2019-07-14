#Figure 2::  Mutual Exclusivity

adeno_AA=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Both_adeno_AA/all_lesions.conf_90.txt', sep='\t')
adeno_EA=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Both_adeno_EA/all_lesions.conf_90.txt', sep='\t')
sq_AA=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Both_sq_AA/all_lesions.conf_90.txt', sep='\t')
sq_EA=read.csv('/home/sinhas8/Projects/Project_Chromotrypsis/4.Results/GISTIC_results/Both_sq_EA/all_lesions.conf_90.txt', sep='\t')

###*HYP
oncoprint_format<-function(Event_Type='Del', input_mat=sq_AA, freq_thresh=0.3, cancer_type='sq', race){
	AA = input_mat[grep(Event_Type, input_mat$Unique.Name),]
	AA = AA[grep('0:',AA$Amplitude.Threshold),]
	df = AA[,c(2, 10:(ncol(AA)-1))]

	df=df[!duplicated(df[,1]),]
	rownames(df)= sapply(df[,1], function(x) paste(gsub(' ','',x), '|', sep=''))
	df=df[,-1]

	df=df[apply(df, 1, function(x) (sum(x>0)/length(x))>freq_thresh ),]
	
	if(Event_Type=='Amp'){shallow='GAIN';Deep='AMP'}
	if(Event_Type=='Del'){shallow='HETLOSS';Deep='HOMDEL'}

	df.list <- split(df, seq(nrow(df)))
	df.list <- setNames(split(df, seq(nrow(df))), rownames(df))
	correct_format_df=t(do.call(cbind,df.list))
	correct_format_df[,1]=gsub('1',shallow,correct_format_df[,1])
	correct_format_df[,1]=gsub('2',Deep,correct_format_df[,1])
	final_df=correct_format_df[correct_format_df[,1]!= '0',]
	to_write=data.frame(t(sapply(names(final_df), function(x) unlist(strsplit(x, '\\|')[[1]])[c(2,1)])), final_df)
	to_write$Type='CNA'
	colnames(to_write)=c('Sample','Gene','Alteration', 'Type')

#write.table(to_write,'/home/sinhas8/Projects/Project_Chromotrypsis/Results_New/adeno_Amp_OncoPrint.csv', sep='\t',quote=F, row.names=F)
	write.table(to_write, paste('/home/sinhas8/Projects/Project_Chromotrypsis/Results_New/',race,'_',cancer_type, Event_Type,  '_OncoPrint.csv', sep=''), sep='\t',quote=F, row.names=F)
}

oncoprint_format(Event_Type='Del', input_mat=adeno_EA, freq_thresh=0.3, cancer_type='ad', race='EA')
oncoprint_format(Event_Type='Amp', input_mat=adeno_EA, freq_thresh=0.3, cancer_type='ad', race='EA')
oncoprint_format(Event_Type='Amp', input_mat=sq_EA, freq_thresh=0.3, cancer_type='sq', race='EA')
oncoprint_format(Event_Type='Del', input_mat=sq_EA, freq_thresh=0.3, cancer_type='sq', race='EA')

