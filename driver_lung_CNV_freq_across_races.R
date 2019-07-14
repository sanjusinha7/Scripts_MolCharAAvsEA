drivers=read.csv()
data.frame( data.frame(drivers$Gene[grep('LUSC',drivers$Cancer)], Type='LUSC'),
            data.frame(drivers$Gene[grep('LUAD',drivers$Cancer)], Type='LUAD'),
            data.frame(drivers$Gene[grep('PANCAN',drivers$Cancer)], Type='PANCAN')
            )
