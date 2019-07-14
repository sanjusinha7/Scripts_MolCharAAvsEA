##For Running GISTIC

cd /home/sinhas/Module/GISTIC
##Format of the command
#./gistic2 -s *Segments_file_address* -refgene *reference_genome* -b *folder to save results* -genegistic 1

##Running all the samples
./gistic2 
-seg /home/sinhas8/Downloads/Segments_OncoScan_try1_infstonegative5_noLOH.txt 
-refgene /home/sinhas8/Modules/GISTIC/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat 
-b /home/sinhas8/GISTIC_Try2_OncoScan_setting_Broad_prostatepaper_genelevel_conf75/  -genegistic 1


##Running just teh Tumor samples
./gistic2 -seg /home/sinhas8/Downloads/Segments_OncoScan_try1_infstonegative5_noLOH_OnlyTumor.txt -refgene /home/sinhas8/Modules/GISTIC/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -b /home/sinhas8/GISTIC_Try2_OncoScan_setting_Broad_prostatepaper_genelevel_conf75_OnlyTumors/  -genegistic 1



###
./gistic2 -seg /home/sinhas8/Downloads/Segments_wdout_LOH_sex_try1_17th_complete_AA.txt -refgene /home/sinhas8/Modules/GISTIC/refgenefiles/hg19.mat -b /home/sinhas8/GISTIC_delthis_AA/  -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1


./gistic2 -seg /home/sinhas8/Downloads/Segments_wdout_LOH_sex_try1_17th_complete_EA.txt -refgene /home/sinhas8/Modules/GISTIC/refgenefiles/hg19.mat -b /home/sinhas8/GISTIC_delthis_EA/  -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1

./gistic2 -seg /home/sinhas8/Downloads/Both_AA.txt -refgene /home/sinhas8/Modules/GISTIC/refgenefiles/hg19.mat -b /home/sinhas8/Both_AA/  -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -twoside 1

./gistic2 -seg /home/sinhas8/Downloads/Both_EA.txt -refgene /home/sinhas8/Modules/GISTIC/refgenefiles/hg19.mat -b /home/sinhas8/Both_EA/  -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -twoside 1



