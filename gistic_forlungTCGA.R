#####Calling GISTIC function
module load gistic
refgenefile=$GISTIC_HOME/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat

###for LUSC AA
mkdir $PWD/LUSC_AA
basedir=$PWD/LUSC_AA
segfile=/home/sinhas8/LUSC_AA.txt
gistic2  -seg $segfile -refgene $refgenefile -b $basedir -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1


###for LUSC EA
mkdir $PWD/LUSC_EA
basedir=$PWD/LUSC_EA
segfile=/home/sinhas8/LUSC_EA.txt
gistic2  -seg $segfile -refgene $refgenefile -b $basedir -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1

###for LUAD AA
mkdir $PWD/LUAD_AA
basedir=$PWD/LUAD_AA
segfile=/home/sinhas8/LUAD_AA.txt
gistic2  -seg $segfile -refgene $refgenefile -b $basedir -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1


###for LUAD EA
mkdir $PWD/LUAD_EA
basedir=$PWD/LUAD_EA
segfile=/home/sinhas8/LUAD_EA.txt
gistic2  -seg $segfile -refgene $refgenefile -b $basedir -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1
get -r LUSC_AA
get -r LUSC_EA
get -r LUAD_AA
get -r LUAD_EA
