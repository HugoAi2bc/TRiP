#! /bin/bash

usage() { echo "Usage: $0 -N <Sample name>" 1>&2 ; echo "\n -N\tsample name\n" ; exit 1; }

while getopts ":N:" option; do
    case "${option}" in
        N)
            N=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if [ -z "${N}" ]; then
    usage
fi

path="/data/RESULTS/qualitativeAnalysis/"
# directory=`ls ${path}`
# for direct in ${directory}
# do
	# cd ${path}/
	# cd ${path}/${direct}

	########### Kmer repartition ###########

# file="${N}.kmerRepartition"
# for file in `ls ${path}kmerRepartition | grep kmerRepartition`
# do
# #Nom de l'échantillon
# sample=${N}
# sample="${file%.*.*}"
echo "#! bin/R" > ${path}graphes/kmerRepartition/${N}.tempoR.R
echo "kmer<-read.table(file = '"${path}"kmerRepartition/"${N}".kmerRepartition.txt')" >> ${path}graphes/kmerRepartition/${N}.tempoR.R
echo "jpeg(filename = '"${path}"graphes/kmerRepartition/"${N}".kmerRepartition.jpeg')" >> ${path}graphes/kmerRepartition/${N}.tempoR.R
echo "barplot(kmer\$V2,names.arg = kmer\$V1)" >> ${path}graphes/kmerRepartition/${N}.tempoR.R
echo "dev.off()" >> ${path}graphes/kmerRepartition/${N}.tempoR.R
R CMD BATCH ${path}graphes/kmerRepartition/${N}.tempoR.R
rm -f ${path}graphes/kmerRepartition/${N}.tempoR.R
# done
# done
