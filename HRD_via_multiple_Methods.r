#########################################################################
#####################     SUPPLEMENTARY METHODS     #####################
#########     Implementation of NtAI, LST and HRD-LOH in R     ##########
#########################################################################


#########################################################################
############################     HOW TO USE    ##########################
#########################################################################

# 'x' is a matrix of segmented output from ASCAT, with at least the
#   following columns (column names are not important):
# 1: sample id
# 2: chromosome (numeric)
# 3: segment start
# 4: segment end
# 5: number of probes
# 6: total copy number
# 7: nA
# 8: nB
# 9: ploidy
# 10: contamination, aberrant cell fraction

chrominfo <- GetChrominfo() # hg19

# chrominfo is a 5 column matrix with information about the chromosomes:
# 1: chromosome number (23 = X, 24 = Y)
# 2: chromosome length
# 3: centromere location (mean of start and end)
# 4: centromere start
# 5: centromere end

ntai <- calc.ai(x, chrominfo)
lst <- calc.lst(x, chrominfo)
hrd <- calc.hrd(x)      # Throughout this text, "hrd" refers to "HRD-LOH"



#########################################################################
############################     FUNCTIONS     ##########################
#########################################################################


#########################################################################
# Number of telomeric AI
# updated and improved to also account for ploidy - 2013-03-11
#########################################################################

calc.ai <- function(
  seg,                # Matrix of segmented output from ASCAT
  chrominfo = FALSE,  # Matrix with information about the chromosomes
  min.size = 0,       # Minimum size of segments
  min.probes = 500,   # Minimum number of probes in segments
  cont = 0,           # Contamination threshold. 0 ignores contamination
  check.names = TRUE, # Rename any duplicated samples
  ploidyByChromosome = TRUE,
  return.loc = FALSE, # Return location of NtAI's
  shrink = TRUE       # Joins segments of identical allelic copy number
  ) {
  if(ploidyByChromosome){cat("Determining chromosome-specific ploidy by major copy number fraction\n")}
  if(!ploidyByChromosome){cat("Determining sample ploidy by major copy number fraction over-all\n")}
  seg[,1] <- as.character(seg[,1])
  if(check.names){
    seg <- check.names.fn(seg)
  }
  # check and potentially fix if nB is always the smaller column
  if(!all(seg[,8] <= seg[,7]) ){
    # In case ASCAT people change the algorithm
    cat("Warning!! nB  not always <= nA!!  -- Correcting for internal use (only!)\n") 
    tmp <- seg
    seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
    seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
  }
  # remove segments smaller min.size and min.probes,
  #   and with too much contamination
  seg <- seg[seg[,5] >= min.probes,]
  seg <- seg[seg[,4]- seg[,3] >= min.size,]
  seg <- seg[seg[,10] >= cont,]
  samples <- as.character(unique(seg[,1]))
  if(shrink){
    new.seg <- seg[1,]
    for(j in samples){
      sample.seg <- seg[seg[,1] %in% j,]
      new.sample.seg <- seg[1,]
      for(i in unique(sample.seg[,2])){
        sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
        sample.chrom.seg <- shrink.seg.ai(sample.chrom.seg)
        new.sample.seg <- rbind(new.sample.seg, sample.chrom.seg)
      }
      new.seg <- rbind(new.seg, new.sample.seg[-1,])    
    }
    seg <- new.seg[-1,]
  }      
  AI <- rep(NA, nrow(seg)) # Add a column to call AI
  seg <- cbind(seg, AI)
  samples <- as.character(unique(seg[,1]))
  ascat.ploidy <- setNames(seg[!duplicated(seg[,1]),9], seg[!duplicated(seg[,1]),1])
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    if(!ploidyByChromosome){
      ploidy <- vector()
      for(k in unique(sample.seg[,6])){
        tmp <- sample.seg[sample.seg[,6] %in% k,]
        ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
      }
      ploidy <- as.numeric(names(ploidy[order(ploidy,decreasing=T)]))[1]    
      # Update "ploidy" column, so the new calculated value can be returned
      sample.seg[,9] <- ploidy 
      # Add a columnm to define AI,
      # with codes for telomeric/interstitial/whole chromosome:
      # 0 = no AI
      # 1 = telomeric
      # 2 = interstitial
      # 3 = whole chromosome
      if(ploidy %in% c(1,seq(2, 200,by=2))){
        sample.seg[,'AI'] <- c(0,2)[match(sample.seg[,7] == sample.seg[,8], c('TRUE', 'FALSE'))]
      }    
      if(!ploidy %in%  c(1,seq(2, 200,by=2))){
        sample.seg[,'AI'] <- c(0,2)[match(sample.seg[,7] + sample.seg[,8] == ploidy & sample.seg[,7] != ploidy, c('TRUE', 'FALSE'))]
      }
    }
    new.sample.seg<- sample.seg[1,]
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(nrow(sample.chrom.seg) == 0){ next}
      if(ploidyByChromosome){
        ploidy <- vector()
        for(k in unique(sample.seg[,6])){
          tmp <- sample.chrom.seg[sample.chrom.seg[,6] %in% k,]
          ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
          ploidy <- ploidy[!names(ploidy) %in% 0] #Remove any ploidy calls of zero
        }
        ploidy <- as.numeric(names(ploidy[order(ploidy,decreasing=T)]))[1]    
        sample.chrom.seg[,9] <- ploidy # update "ploidy" column, so the new calculated value can be returned
        if(ploidy %in% c(1,seq(2, 200,by=2))){
          sample.chrom.seg[,'AI'] <- c(0,2)[match(sample.chrom.seg[,7] == sample.chrom.seg[,8], c('TRUE', 'FALSE'))]
        }    
        if(!ploidy %in%  c(1,seq(2, 200,by=2))){
          sample.chrom.seg[,'AI'] <- c(0,2)[match(sample.chrom.seg[,7] + sample.chrom.seg[,8] == ploidy & sample.chrom.seg[,8] != 0, c('TRUE', 'FALSE'))]
        }
        sample.seg[sample.seg[,2] %in% i,9] <-ploidy
        sample.seg[sample.seg[,2] %in% i,'AI'] <-sample.chrom.seg[,'AI']
      }
      # By logical, we assume that chrominfo == FALSE, hence we here proceed without checking for the centromere (useful for non-human samples)
      if(class(chrominfo) == 'logical'){      
        if(sample.chrom.seg[1,'AI'] == 2 & nrow(sample.chrom.seg) != 1){
          sample.seg[sample.seg[,2]==i,'AI'][1] <- 1
        }
        if(sample.chrom.seg[nrow(sample.chrom.seg),'AI'] == 2 & nrow(sample.chrom.seg) != 1){
          sample.seg[sample.seg[,2]==i,'AI'][nrow(sample.seg[sample.seg[,2]==i,])] <- 1
        }
      }
      if(class(chrominfo) != 'logical'){    # Here we consider the centromere    
        if(sample.chrom.seg[1,'AI'] == 2 & nrow(sample.chrom.seg) != 1 & sample.chrom.seg[1,4] < (chrominfo[i,3])){
          sample.seg[sample.seg[,2]==i,'AI'][1] <- 1
        }
        if(sample.chrom.seg[nrow(sample.chrom.seg),'AI'] == 2 & nrow(sample.chrom.seg) != 1 & sample.chrom.seg[nrow(sample.chrom.seg),3] > (chrominfo[i,3])){
          sample.seg[sample.seg[,2]==i,'AI'][nrow(sample.seg[sample.seg[,2]==i,])] <- 1
        }
      }
      if(nrow(sample.seg[sample.seg[,2]==i,]) == 1 & sample.seg[sample.seg[,2]==i,'AI'][1] != 0){
        sample.seg[sample.seg[,2]==i,'AI'][1] <- 3
      }
    }
    seg[seg[,1] %in% j,] <- sample.seg
  }  
  samples <- as.character(unique(seg[,1]))
  no.events <- matrix(0, nrow=length(samples), ncol=12)
  rownames(no.events) <- samples
  colnames(no.events) <- c("Telomeric AI", "Mean size", "Interstitial AI", "Mean Size", "Whole chr AI", "Telomeric LOH",  "Mean size", "Interstitial LOH", "Mean Size", "Whole chr LOH", "Ploidy", "Aberrant cell fraction")
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    no.events[j,1] <- nrow(sample.seg[sample.seg[,'AI'] == 1,])
    no.events[j,2] <- mean(sample.seg[sample.seg[,'AI'] == 1,4] - sample.seg[sample.seg[,'AI'] == 1,3])
    no.events[j,3] <- nrow(sample.seg[sample.seg[,'AI'] == 2,])
    no.events[j,4] <- mean(sample.seg[sample.seg[,'AI'] == 2,4] - sample.seg[sample.seg[,'AI'] == 2,3])
    no.events[j,5] <- nrow(sample.seg[sample.seg[,'AI'] == 3,])
    no.events[j,11] <- ascat.ploidy[j]
    no.events[j,12] <- unique(sample.seg[,10]) # aberrant cell fraction
    # Here we restrict ourselves to real LOH
    sample.seg <- sample.seg[sample.seg[,8] == 0,]
    no.events[j,6] <- nrow(sample.seg[sample.seg[,'AI'] == 1,])
    no.events[j,7] <- mean(sample.seg[sample.seg[,'AI'] == 1,4] - sample.seg[sample.seg[,'AI'] == 1,3])
    no.events[j,8] <- nrow(sample.seg[sample.seg[,'AI'] == 2,])
    no.events[j,9] <- mean(sample.seg[sample.seg[,'AI'] == 2,4] - sample.seg[sample.seg[,'AI'] == 2,3])
    no.events[j,10] <- nrow(sample.seg[sample.seg[,'AI'] == 3,])
  }
  if(return.loc){
    return(seg)
  } else {
    return(no.events)
  }
}




#########################################################################
# This function calculates Popova's LST measure
# (Popova 2012, Cancer research)
# http://www.bio-protocol.org/e814
# Popova's cutoffs: 15 LSTs in near-diploid, 20 in near-tetraploid
#########################################################################

calc.lst <- function(
  seg,                # Matrix of segmented output from ASCAT
  chrominfo,          # Matrix with information about the chromosomes
  nA = 7,             # The column index of 'seg' where copy number of A
                      #   allele is found
  check.names = T,    # Rename any duplicated samples
  min.probes = 50,    # As described by Popova (50 for SNP6)  
  return.loc = FALSE, # Return location of LST sites
  chr.arm = 'no'      # Option to use chromosome arms defined during
                      #   segmentation. The option must give a column
                      #   that holds the chromosome arm information
                      #   (or 'no' to not use this)
  ){
  seg[,1] <- as.character(seg[,1])
  if(check.names){
    seg <- check.names.fn(seg)
  }
  if(! all(seg[,8] <= seg[,7]) ){
    # In case ASCAT people change the algorithm  ### They did!!
    cat("Warning!! nB  not always <= nA!!  -- Correcting for internal use (only!)\n") 
    tmp <- seg
    seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
    seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
  }
  seg <- seg[seg[,5] >= min.probes,]
  nB <- nA+1
  samples <- unique(seg[,1])
  output <- setNames(rep(0,length(samples)), samples)
  if(return.loc) {
      out.seg <- matrix(0,0,10)
      colnames(out.seg) <- c(colnames(seg)[1:8],'LST breakpoint', 'chr. arm')
  }
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    sample.lst <- c()
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(chr.arm !='no'){
        p.max <- if(any(sample.chrom.seg[,chr.arm] == 'p')){max(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'p',4])}
        q.min <- min(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'q',3])
      }
      if(nrow(sample.chrom.seg) < 2) {next}
      sample.chrom.seg.new <- sample.chrom.seg
      if(chr.arm == 'no'){
      	# split into chromosome arms
        p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,3] <= chrominfo[i,4],] 
        q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,4] >= chrominfo[i,5],]
        q.arm<- shrink.seg.ai(q.arm)
        p.arm<- shrink.seg.ai(p.arm)
        p.arm[nrow(p.arm),4] <- chrominfo[i,4]
        q.arm[1,3] <- chrominfo[i,5]
      }
      if(chr.arm != 'no'){
        q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'q',,drop=F]
        q.arm<- shrink.seg.ai(q.arm)
        q.arm[1,3] <- q.min
        if(any(sample.chrom.seg.new[,chr.arm] == 'p')){
          # split into chromosome arms
          p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'p',,drop=F] 
          p.arm<- shrink.seg.ai(p.arm)
          p.arm[nrow(p.arm),4] <- p.max
        }
      }  
      n.3mb <- which((p.arm[,4] - p.arm[,3]) < 3e6)
      while(length(n.3mb) > 0){
        p.arm <- p.arm[-(n.3mb[1]),]
        p.arm <- shrink.seg.ai(p.arm)
        n.3mb <- which((p.arm[,4] - p.arm[,3]) < 3e6)
      }
      if(nrow(p.arm) >= 2){
        p.arm <- cbind(p.arm[,1:8], c(0,1)[match((p.arm[,4]-p.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
        for(k in 2:nrow(p.arm)){
          if(p.arm[k,9] == 1 & p.arm[(k-1),9]==1 & (p.arm[k,3]-p.arm[(k-1),4]) < 3e6){
            sample.lst <- c(sample.lst, 1)
            if(return.loc){
              ## Number indicates if start (1) or end (2) defines the breakpoint
              a<- cbind(p.arm[(k-1),1:8], 2,'p-arm') 
              b <- cbind(p.arm[k,1:8], 1,'p-arm')
              colnames(a)[9:10]<- colnames(b)[9:10]<- c('LST breakpoint', 'chr. arm')
              out.seg <- rbind(out.seg, a,b)
            }
          }
        }
      }
      n.3mb <- which((q.arm[,4] - q.arm[,3]) < 3e6)
      while(length(n.3mb) > 0){
        q.arm <- q.arm[-(n.3mb[1]),]
        q.arm <- shrink.seg.ai(q.arm)
        n.3mb <- which((q.arm[,4] - q.arm[,3]) < 3e6)
      }
      if(nrow(q.arm) >= 2){
        q.arm <- cbind(q.arm[,1:8], c(0,1)[match((q.arm[,4]-q.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
        for(k in 2:nrow(q.arm)){
          if(q.arm[k,9] == 1 & q.arm[(k-1),9]==1 & (q.arm[k,3]-q.arm[(k-1),4]) < 3e6){
            sample.lst <- c(sample.lst, 1)
            if(return.loc){
              ## Number indicates if start (1) or end (2) defines the breakpoint
              a<- cbind(q.arm[(k-1),1:8], 2,'q-arm') 
              b <- cbind(q.arm[k,1:8], 1,'q-arm')
              colnames(a)[9:10]<- colnames(b)[9:10]<- c('LST breakpoint', 'chr. arm')
              out.seg <- rbind(out.seg, a,b)
            }
          }
        }
      }
    }
    output[j] <- sum(sample.lst)
  }
  if(return.loc){
    return(out.seg)
  } else {
    return(output)
  }
}
    


#########################################################################
# This function is an implementation of the cisplatin predictor developed
# by Myriad with Gordon Mills PMID: 23047548.
# Abkevichs cutoffs: > 10, found in the supplementary info
#########################################################################

calc.hrd <- function(
  seg,                # Matrix of segmented output from ASCAT
  nA = 7,             # The column index of 'seg' where copy number of A
                      #   allele is found
  check.names = TRUE, # Rename any duplicated samples
  return.loc = FALSE  # Return location of HRD sites
  ){
  seg[,1] <- as.character(seg[,1])
  if(check.names){
    seg <- check.names.fn(seg)
  }
  if(! all(seg[,8] <= seg[,7]) ){
    # In case ASCAT people change the algorithm  ### They did!!
    cat("Warning!! nB  not always <= nA!!  -- Correcting for internal use (only!)\n") 
    tmp <- seg
    seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
    seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
  }
  nB <- nA+1
  output <- rep(0, length(unique(seg[,1])))
  names(output) <- unique(seg[,1])
  if(return.loc) {
      out.seg <- matrix(0,0,9)
      colnames(out.seg) <- c(colnames(seg)[1:8],'HRD breakpoint')
  }
  for(i in unique(seg[,1])){
    segSamp <- seg[seg[,1] %in% i,]
    chrDel <-vector() 
    for(j in unique(segSamp[,2])){   
      if(all(segSamp[segSamp[,2] == j,nB] == 0)) {
        chrDel <- c(chrDel, j)
      }
    }
    segLOH <- segSamp[segSamp[,nB] == 0 & segSamp[,nA] != 0,,drop=F]
    segLOH <- segLOH[segLOH[,4]-segLOH[,3] > 15e6,,drop=F]
    segLOH <- segLOH[!segLOH[,2] %in% chrDel,,drop=F]
    output[i] <- nrow(segLOH)
    if(return.loc){
      if(nrow(segLOH) < 1){next}
      segLOH <- cbind(segLOH[,1:8], 1)
      colnames(segLOH)[9] <- 'HRD breakpoint'
      out.seg <- rbind(out.seg, segLOH)
    }
  }
  if(return.loc){
    return(out.seg)
  } else {
    return(output)
  }
}  



#########################################################################
# Function to check for duplicate sample names in an ASCAT segmented
# object, and return new names
#########################################################################

check.names.fn <- function(
  seg,
  max.size = 3.2e9,
  remove.dup = TRUE
  ){
  tmp <- setNames(paste(seg[,1],seg[,9],seg[,11],sep='_'), seg[,1])
  sample.names <-  names(tmp[!duplicated(tmp)]) # this is based on the fact that ASCAT ploidy is always highly variable, as it is a mean of all the segments
  if(any(duplicated(sample.names))){
    cat("Warning!! Duplicate names found! Renaming duplicates, adding a '.(number)' to each beyond no. 1\n")
    dup.samp <- unique(sample.names[duplicated(sample.names)])
    for(i in 1:length(dup.samp)){
      dup.samp.seg <- seg[seg[,1] %in% dup.samp[i],]    
      tmp<- unique(paste(dup.samp.seg[,1],dup.samp.seg[,9],dup.samp.seg[,11],sep='_'))
      for(k in 2:length(tmp)){
        dup.samp.seg[paste(dup.samp.seg[,1],dup.samp.seg[,9],dup.samp.seg[,11],sep='_') %in% tmp[k],1] <- paste(dup.samp.seg[paste(dup.samp.seg[,1],dup.samp.seg[,9],dup.samp.seg[,11],sep='_') %in% tmp[k],1],'.',k,sep='')
      }
      seg[seg[,1] %in% dup.samp[i],] <- dup.samp.seg
    }
  }
  # Now check for "silent" duplicates
  sample.names <- unique(seg[,1])
  for(i in sample.names){
    tmp <- seg[seg[,1] %in% i,]
    if(sum(tmp[,4] - tmp[,3]) > max.size){
      a <- which(tmp[,2] %in% 1)
      b <- which(tmp[,2] %in% 22)
      a <- a[!(a == (1:b[1])[1:length(a)])][1]
      tmp[a:nrow(tmp),1] <- paste('Duplicate', tmp[1,1], sep='-')
      seg[seg[,1] %in% i,] <- tmp
    }
  }
  if(remove.dup){
    seg <- seg[!grepl('Duplicate', seg[,1]),]
  }
  return(seg)
}



#########################################################################
# AI-specific function to shrink a DNAcopy style matrix.
# Works on a chromosome level.
# Used to condense rows with identical CN values, which may have been
# separated due to other values not currently of interest (e.g. filtered
# out for number of probes)
#########################################################################

shrink.seg.ai <- function(chr.seg)  {  
  new.chr <- matrix(0,0,ncol(chr.seg))
  colnames(new.chr) <- colnames(chr.seg)
  new.chr <- chr.seg
  seg.class <- c(1)
  for(j in 2:nrow(new.chr)){
    ifelse(new.chr[(j-1),7] == new.chr[j,7] & new.chr[(j-1),8] == new.chr[j,8], seg.class <- c(seg.class, seg.class[j-1]), seg.class <- c(seg.class, seg.class[j-1]+1))
  }
  for(j in unique(seg.class)){
    new.chr[seg.class %in% j,4] <- max(new.chr[seg.class %in% j,4])
    new.chr[seg.class %in% j,5] <- sum(new.chr[seg.class %in% j,5])
  }
  new.chr<- new.chr[!duplicated(seg.class),]
  return(new.chr)
}

#########################################################################
# Functions used to get/make the chrominfo table, which is needed to
# calculate NtAI and LST scores.
#########################################################################

GetGzFromUrl <- function(url, ...) {
  # http://stackoverflow.com/questions/9548630/read-gzipped-csv-directly-from-a-url-in-r
  con <- gzcon(url(url))
  txt <- readLines(con)
  dat <- read.delim(textConnection(txt), ...)  
  return(dat)
}

GetChrominfo <- function() {
  # Get chromInfo table from UCSC
  chrom <- GetGzFromUrl("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz", header = FALSE)
  chrom <- subset(chrom, grepl("^chr[0-9XY]{1,2}$", chrom[,1]))
  # Get gap table from UCSC
  gaps <- GetGzFromUrl("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz", header = FALSE)
  centro <- subset(gaps, gaps[,8] == "centromere")
  # Merge the relevant info from the two tables
  chrominfo <- merge(chrom[,1:2], centro[,2:4], by.x = 1, by.y = 1) # merge size and centromere location
  chrominfo$centromere <- rowMeans(chrominfo[,3:4]) # convert centromere start and end into one location (the mean)
  chrominfo <- chrominfo[,c(1,2,5,3,4)] # keep chromosome, size and centromere location
  colnames(chrominfo) <- c("chr", "size", "centromere", "centstart", "centend")
  chrominfo[,1] <- as.character(chrominfo[,1])
  chrominfo$chr <- sub("chr", "", chrominfo$chr)
  chrominfo$chr <- sub("X", "23", chrominfo$chr)
  chrominfo$chr <- sub("Y", "24", chrominfo$chr)
  chrominfo[,1] <- as.numeric(chrominfo[,1])
  chrominfo <- chrominfo[order(chrominfo$chr), ]  # order by chromosome number
  rownames(chrominfo) <- as.character(chrominfo[,1])
  chrominfo <- as.matrix(chrominfo)
  return(chrominfo)  
}
