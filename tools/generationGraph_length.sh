#! /bin/bash

path="/data/RESULTS/qualitativeAnalysis"
# directory=`ls ${path}`
# for direct in ${directory}
# do
	# cd ${path}/
	# cd ${path}/${direct}
	
	########### Kmer repartition ###########

	for file in `ls ${path}/kmerRepartition | grep kmerRepartition`
	do
		#Nom de l'échantillon
		sample="${file%.*.*}"
		echo "#! bin/R" > ${path}/graphes/kmerRepartition/tempoR.R
		echo "kmer<-read.table(file = '"${path}"/kmerRepartition/"${file}"')" >> ${path}/graphes/kmerRepartition/tempoR.R
		echo "jpeg(filename = '"${path}"/graphes/kmerRepartition/"${sample}.kmerRepartition.jpeg"')" >> ${path}/graphes/kmerRepartition/tempoR.R
		echo "barplot(kmer\$V2,names.arg = kmer\$V1)" >> ${path}/graphes/kmerRepartition/tempoR.R
		echo "dev.off()" >> ${path}/graphes/kmerRepartition/tempoR.R
		R CMD BATCH ${path}/graphes/kmerRepartition/tempoR.R
		rm -f ${path}/graphes/kmerRepartition/tempoR.R
	done
# done
