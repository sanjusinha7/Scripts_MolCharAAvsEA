####Genotyping
require(crlmm)
BiocManager::install('hapmapsnp6')
path <- system.file("celFiles", package="hapmapsnp6")
celFiles <- list.celfiles('/Users/sinhas8/Downloads/CS_022389_BridRyan_run7_CEL/', full.names=TRUE)
celFiles=celFiles[grep('C\\.',celFiles)]
system.time(crlmmResult <- crlmm(celFiles, verbose=TRUE))
