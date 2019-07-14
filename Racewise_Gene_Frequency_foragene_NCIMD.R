##Below are frequency files.
err_handle<-function(x){ tryCatch(x, error=function(e){NA}) }
setwd('/Users/sinhas8/')
Freq_LUSC=read.csv('./Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Freq_LUSC.csv')
Freq_LUAD=read.csv('./Project_Chromotrypsis/4.Results/GISTIC_results/Trimmed_GISTIC_Results/Freq_LUAD.csv')
###Starting with ~3400 sq and ~2500 candidate recurrent genes
dim(Freq_LUSC);dim(Freq_LUAD)
Candidate_AASpecific_Driver<-function(df=Freq_LUAD, cancer_type='LUAD', K=1.5){
  df=df[!(df$Aberrance_AA !='' &df$Aberrance_EA !=''),]
  df=df[((df$Aberrance_AA=='Del' & df$Del_AA_freq>(K*df$Del_EA_freq)) & 
           (df$Del_EA_freq+df$Del_AA_freq)>0.10 )|
          ((df$Aberrance_AA=='Amp' & df$Amp_AA_freq>(K*df$Amp_EA_freq)) & 
             (df$Del_EA_freq+df$Del_AA_freq)>0.10),]
  df
}

Freq_LUSC=data.frame(Candidate_AASpecific_Driver(df=Freq_LUSC, cancer_type='LUSC', K=1.5), hist='LUSC')
Freq_LUAD=data.frame(Candidate_AASpecific_Driver(df=Freq_LUAD, cancer_type='LUAD', K=1.5), hist='LUAD')
dim(Freq_LUSC); dim(Freq_LUAD)
###############
##Adding Survival
###############
CNV_FORSurvival=CNV[,na.omit(match(mat$acc[mat$hist=='LUAD'], colnames(CNV)))]
mat_FORSurvival=mat[!is.na(match(mat$acc, colnames(CNV_FORSurvival))),]
temp=apply(CNV_FORSurvival[Driver_OI,], 1, function(x) 
  summary(coxph(Surv(mat_FORSurvival$survival, mat_FORSurvival$lungcancer_death_5_years)~ x))$coefficients[c(1,5)])
