#! /bin/bash

usage() { echo "Usage: $0 -N <Sample name> -l <Read length> -m <before start position> -M <after start position>" 1>&2 ; echo "\n -N\tsample name\n -l\tRead length\n -m <int>\tnumber of base before start codon\n -M <int>\tnumber of base after start codon\n" ; exit 1; }

while getopts ":N:l:m:M:" option; do
    case "${option}" in
        N)
            N=${OPTARG}
            ;;
        l)
            l=${OPTARG}
            ;;
        m)
            m=${OPTARG}
            ;;
        M)
			M=${OPTARG}
			;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if [ -z "${N}" ] || [ -z "${l}" ] || [ -z "${m}" ] || [ -z "${M}" ]; then
    usage
fi

path="/data/RESULTS/qualitativeAnalysis/"

########### Periodicity ###########
for pos in start stop
do
	if [ ${pos} = "start" ]; then
		file="${N}.${l}.periodicity.start.CDS.-${m}+${M}"
	else
		file="${N}.${l}.periodicity.stop.CDS.-${M}+${m}"
	fi

	echo "#! bin/R" > "${path}graphes/periodicity/${N}.${l}.${m}.${M}.tempoR.R"
	echo "perio<-read.table(file = '${path}periodicity/${file}.txt')" >> "${path}graphes/periodicity/${N}.${l}.${m}.${M}.tempoR.R"
	echo "jpeg(filename = '${path}graphes/periodicity/${file}.jpeg')" >> "${path}graphes/periodicity/${N}.${l}.${m}.${M}.tempoR.R"
	echo "barplot(perio\$V2,col = c('red','green','blue'),names.arg = perio\$V1,cex.names = 0.75, las=3, xlab = 'Relative position of reads start', ylab = 'Number of reads')" >> "${path}graphes/periodicity/${N}.${l}.${m}.${M}.tempoR.R"
	echo "dev.off()" >> "${path}graphes/periodicity/${N}.${l}.${m}.${M}.tempoR.R"
	R CMD BATCH "${path}graphes/periodicity/${N}.${l}.${m}.${M}.tempoR.R"
	rm -f "${path}graphes/periodicity/${N}.${l}.${m}.${M}.tempoR.R"
done
