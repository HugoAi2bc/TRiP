#! /bin/bash

path="/home/hugo.arbes/Hugo/GST/ChemoR/prerun/results/qualitativeAnalysis/"
cd ${path}

########### Kmer repartition ###########
mkdir -p ./graphes/kmerRepartition

for file in `ls kmerRepartition | grep kmerRepartition`
do
	#Nom de l'échantillon
	sample="${file%.*.*}"

	echo $sample

	echo "#! bin/R" > tempoR.R
	echo "kmer<-read.table(file = '"${path}"kmerRepartition/"${file}"')" >> tempoR.R
	echo "jpeg(filename = '"${path}graphes/kmerRepartition/${sample}.kmerRepartition.jpeg"')" >> tempoR.R
	echo "barplot(kmer\$V2,names.arg = kmer\$V1)" >> tempoR.R
	echo "dev.off()" >> tempoR.R
	R CMD BATCH tempoR.R
	rm -f tempoR.R
done

# ########## Periodicity ###########
mkdir -p ./graphes/periodicity

for file in `ls periodicity`
do
	#Nom de l'échantillon
	sample="${file%.*}"
	echo $sample
	echo "#! bin/R" > tempoR.R
	echo "perio<-read.table(file = '"${path}"periodicity/"${file}"')" >> tempoR.R
	echo "jpeg(filename = '"${path}"graphes/periodicity/"${sample}".jpeg')" >> tempoR.R
	echo "barplot(perio\$V2,col = c('red','green','blue'),names.arg = perio\$V1,cex.names = 0.75, las=3)" >> tempoR.R
	echo "dev.off()" >> tempoR.R
	R CMD BATCH tempoR.R
	rm -f tempoR.R
done
