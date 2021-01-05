#! /bin/bash

path="/data/RESULTS/qualitativeAnalysis"
# directory=`ls ${path}`

# for direct in ${directory}
# do
	# cd ${path}/
	# cd ${path}/${direct}
	mkdir ${path}/graphes

	########### Periodicity ###########
	mkdir ${path}/graphes/periodicity

	for file in `ls ${path}/periodicity`
	do
		#Nom de l'échantillon
		sample="${file%.*}"
		echo "#! bin/R" > ${path}/graphes/periodicity/tempoR.R
		echo "perio<-read.table(file = '"${path}"/periodicity/"${file}"')" >> ${path}/graphes/periodicity/tempoR.R
		echo "jpeg(filename = '"${path}/"graphes/periodicity/"${sample}".jpeg')" >> ${path}/graphes/periodicity/tempoR.R
		echo "barplot(perio\$V2,col = c('red','green','blue'),names.arg = perio\$V1,cex.names = 0.75, las=3)" >> ${path}/graphes/periodicity/tempoR.R
		echo "dev.off()" >> ${path}/graphes/periodicity/tempoR.R
		R CMD BATCH ${path}/graphes/periodicity/tempoR.R
		rm -f ${path}/graphes/periodicity/tempoR.R
	done
# done
