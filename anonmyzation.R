prev=readxl::read_xlsx('/Users/sinhas8/Downloads/ryan_dataset/Sample ID Legend.xlsx')
new_mapping=readxl::read_xlsx('/Users/sinhas8/Downloads/Mapping_Updated.xlsx', col_names = F)
annotation=read.csv('/Users/sinhas8/Project_Chromotrypsis/2.Data/all_features_withannotation_1oct.csv')

folder = "/Users/sinhas8/Project_Chromotrypsis/2.Data/All_OSCHP_data/"
setwd(folder)
files <- list.files(full.names = F)
files_acc=gsub('[tT]','',sapply(as.character(files), function(x) strsplit(x, '\\.')[[1]][1]))
files_acc[is.na(as.numeric(files_acc))]=sapply(files_acc[is.na(as.numeric(files_acc))], function(x) strsplit(x, '_')[[1]][2])
files_acc[is.na(as.numeric(files_acc))]=sapply(gsub('[A-z]','',names(files_acc[is.na(as.numeric(files_acc))])),
       function(x) strsplit(x, '\\.')[[1]][1])
files_ofInterest=files[match(annotation$acc, files_acc)]


head(files_ofInterest)
# Comparison to the previous study
sum(!is.na(match(annotation$acc, prev$`Patient ID`)))
previous_samples_mapping=prev[!is.na(match(prev$`Patient ID`, annotation$acc)),1:2]
previous_samples_mapping$fileName=names(files_acc[match(previous_samples_mapping$`Patient ID`, files_acc)])
newSamples_df=data.frame(C=paste('Patient', 172:(172+114), sep='_'),
           B=files_acc[match(annotation$acc, files_acc)][is.na(match(files[match(annotation$acc, files_acc)], previous_samples_mapping$fileName))],
           A=files[match(annotation$acc, files_acc)][is.na(match(files[match(annotation$acc, files_acc)], previous_samples_mapping$fileName))])
previous_samples_mapping=data.frame(previous_samples_mapping)
colnames(newSamples_df)=colnames(previous_samples_mapping)
df=rbind(previous_samples_mapping, newSamples_df)
K=1
sapply(1:length(df$fileName),FUN=function(K){
  file.rename(from=df$SRA.ID[K],to= paste(df$SRA.ID[K], '.OSCHP', sep='') )
})
sapply(files[is.na(match(files_acc, annotation$acc))], unlink)
# Samples DS
write.csv(data.frame(df$SRA.ID, paste(df$SRA.ID, '.OSCHP', sep='')), '/Users/sinhas8/Downloads/ryan_dataset/Phenotype_Data_Update/DS_newPart.csv', 
          quote=F)
# Subject DS
annotation$SRA_ID=df$SRA.ID[match(annotation$acc, df$Patient.ID)]
df2writeSample_DS=data.frame(newSamples_df$SRA.ID, ncimd[match(newSamples_df$Patient.ID, ncimd$acc),
                                                         c('smoke', 'age',  'race', 'GENDER','packyrs')])
write.csv(df2writeSample_DS, '/Users/sinhas8/Downloads/ryan_dataset/Phenotype_Data_Update/Sub_DS.csv',
          quote=F)
df2writeSample_Attributes=data.frame(SAMPLE_ID=paste(df$SRA.ID, '.OSCHP', sep=''), 
                                     BODY_SITE='Lung',
                                     ANALYTE_TYPE='DNA',
                                     IS_TUMOR='Y',
                                     HISTOLOGICAL_TYPE=mat$hist[match(df$Patient.ID, mat$acc)],
                                     TUMOR_STAGE=mat$stage[match(df$Patient.ID, mat$acc)],
                                     SEQUENCING_PANEL='OncoVar_CNV')

write.csv(df2writeSample_Attributes, '/Users/sinhas8/Downloads/ryan_dataset/Phenotype_Data_Update/Samples_DS.csv',
          quote=F)


# Updating the mistake
Sub_DS=read.csv('/Users/sinhas8/Downloads/ryan_dataset/Phenotype_Data_Update/5a_SubjectPhenotypes_DS.txt', sep='\t')
levels(Sub_DS$RACE)=c(levels(Sub_DS$RACE)[2],levels(Sub_DS$RACE)[2], 
                      levels(Sub_DS$RACE)[4],levels(Sub_DS$RACE)[4])
levels(Sub_DS$SEX)=c(levels(Sub_DS$SEX)[2],levels(Sub_DS$SEX)[2], 
                     levels(Sub_DS$SEX)[4],levels(Sub_DS$SEX)[4])
write.table(Sub_DS, '/Users/sinhas8/Downloads/ryan_dataset/Phenotype_Data_Update/Subject_DS_Upd.txt', sep='\t', quote = F, row.names = F)
