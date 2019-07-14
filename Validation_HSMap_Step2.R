###After finding modulators for all the pairs which are active
Test_a_GenePair_Modulators_Direction<-function(GenePair, filelist=All_files){

	Mod_Score = readRDS(filelist[ grepl(as.character(unlist(GenePair[1])),filelist) & 
                                     grepl(as.character(unlist(GenePair[2])),filelist)][1])
	Modulators_forthe_Pair = Threshold(Mod_Score, Thresh_EB=0.1, Thresh_MB=0.1, Thresh_CB=0.1)

	if(!is.null(nrow(Modulators_forthe_Pair[[1]])) ){
		Exp_Mods  = Modulators_forthe_Pair[[1]][,c(1,4)]
		Exp_Score = sum(apply(Exp_Mods, 1, function(x) Testing_Modulators_Direction
                               (Mod_GeneName=as.character(x), Type='Exp'))*(nrow(Exp_Mods):1), na.rm=T)
		Exp_here  = c(Exp_Score= Exp_Score, Ex_Exp_Score=(nrow(Exp_Mods)*(nrow(Exp_Mods)+1))/4)
	} else{ Exp_here  = c(Exp_Score=NA, Ex_Exp_Score=NA)}

	if(!is.null(nrow(Modulators_forthe_Pair[[3]])) ){
		CNV_Mods  = Modulators_forthe_Pair[[3]][,c(1,4)]
		CNV_Score = sum(apply(CNV_Mods, 1, Testing_Modulators_Direction)*(nrow(CNV_Mods):1), na.rm=T)
		CNV_here  = c(CNV_Score= CNV_Score, Ex_CNV_Score=(nrow(CNV_Mods)*(nrow(CNV_Mods)+1))/4)
	} else{ CNV_here=c(CNV_Score=NA, Ex_CNV_Score=NA)}
	c(Exp_here, CNV_here)
}

###TOp K Modulators for each pair
Test_a_GenePair_TopK_Modulators_Direction<-function(GenePair, TopK=10, filelist= All_files){

	Mod_Score = readRDS(filelist[ grepl(as.character(unlist(GenePair[1])),filelist) & 
                                     grepl(as.character(unlist(GenePair[2])),filelist)][1])
	Modulators_forthe_Pair = Top_K_Modulators(Mod_Score, K=TopK)

	Exp_Mods  = Modulators_forthe_Pair[[1]][,c(1,4)]
	Exp_Score = sum(apply(Exp_Mods, 1, function(x) Testing_Modulators_Direction
                       (Mod_GeneName=as.character(x), Type='Exp'))*(nrow(Exp_Mods):1), na.rm=T)
	Exp_here  = c(Exp_Score= Exp_Score, Ex_Exp_Score=(nrow(Exp_Mods)*(nrow(Exp_Mods)+1))/4)

	CNV_Mods  = Modulators_forthe_Pair[[3]][,c(1,4)]
	CNV_Score = apply(CNV_Mods, 1, Testing_Modulators_Direction)*(nrow(CNV_Mods):1)
	
	CNV_here  = c( CNV_Score= sum(CNV_Score, na.rm=T),
                       Ex_CNV_Score=(sum(!is.na(CNV_Score))*(sum(!is.na(CNV_Score))+1))/4)

	c(Exp_here, CNV_here)
}

##Given a modulator Gene and a pair with it's score.
Testing_Modulators_Direction<-function(Mod_GeneName, Type='CNV', Percentile='Top'){
	if(Type=='CNV') { test_mat=CCLE_CNV 
			} else{test_mat=CCLE_Exp}

	Expected_direction = ifelse(Mod_GeneName[2]>0, 'Pos', 'Neg')
	ModName=as.character(unlist(Mod_GeneName[1]))
	if(Expected_direction=='Pos'){
		obj_return = tryCatch(test_mat[ModName, grep('K562', colnames(test_mat))] >
                             test_mat[ModName, grep('JURKAT', colnames(test_mat))], error=function(err){NA})
	}	else{
		obj_return = tryCatch(test_mat[ModName, grep('K562', colnames(test_mat))] <
                             test_mat[ModName, grep('JURKAT', colnames(test_mat))], error=function(err){NA})
	}
	obj_return
}

##
