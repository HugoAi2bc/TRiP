#!/bin/sh

usage() { echo "Usage: $0 -m <Shortest kmer length> -M <Longest kmer length> -D <path to BAM by length> -N < Sample name> -F <Fasta file> -O <output dir>" 1>&2 ; echo "\n -N\tsample name\n -D\tpath to directory with uniq and sort BAM files\n -m\tShorter kmer (reads) length\n -M\tLongest kmer (reads) length\n -F\tGenome sequence in fasta format\n -O\tOutput directory" ; exit 1; }

while getopts ":D:m:M:N:F:O:" option; do
    case "${option}" in
        D)
            D=${OPTARG}
            ;;
        m)
            m=${OPTARG}
            ;;
        M)
            M=${OPTARG}
            ;;
        N)
            N=${OPTARG}
            ;;
        F)
            F=${OPTARG}
            ;;
        O)
            O=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if [ -z "${D}" ] || [ -z "${m}" ] || [ -z "${M}" ] || [ -z "${N}" ] || [ -z "${F}" ]; then
    usage
fi
if [ -z "${O}" ]; then
    O="."
fi

mkdir -p ${O}"kmerRepartition"
mkdir -p ${O}"bedCount"
mkdir -p ${O}"sequenceBedCount"

#K-mer repartition
for kmer in `seq ${m} ${M}`;
do
  	echo "Analyse kmer" ${kmer};
  	bamToBed -i ${D}${N}.${kmer}.uniq.sort.bam > ${O}kmerRepartition/${N}.${kmer}.bed ;
  	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' ${O}kmerRepartition/${N}.${kmer}.bed | uniq -c | tr -s ' ' | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$2":"$3"-"$4"("$5")\t"$5}' | sort -k 5 > ${O}bedCount/${N}.${kmer}.count.bed;
  	bedtools getfasta -fi ${F} -bed ${O}bedCount/${N}.${kmer}.count.bed -s -tab -fullHeader -fo ${O}bedCount/${N}.${kmer}.fa;
  	uniq ${O}bedCount/${N}.${kmer}.fa > ${O}bedCount/${N}.${kmer}.uniq.fa;
  	join -1 5 -2 1 -o 1.1,1.2,1.3,1.4,2.2,1.6 ${O}bedCount/${N}.${kmer}.count.bed ${O}bedCount/${N}.${kmer}.uniq.fa | tr " " "\t" > ${O}sequenceBedCount/${N}.${kmer}.count.sequence.bed;
  	rm -f ${O}bedCount/${N}.${kmer}.fa ${O}bedCount/${N}.${kmer}.uniq.fa;
done
