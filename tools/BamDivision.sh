#! /bin/bash

###
###This script allowed the division of SAM file according to kmer length and turns it into BAM
###

usage() { echo "Usage: $0 -m <Shortest kmer length> -M <Longest kmer length> -S <SAM file> -T <Threads> -N <sample name> -O <output dir>" 1>&2 ; echo -e "\n -S\tSAM files (with path) without multimapping and unalign reads\n -m\tShorter kmer (reads) length\n -M\tLongest kmer (reads) length\n -T\tThreads number allows (default value : 4)\n -N\tSample name\n -O\tOutput directory" ; exit 1; }

#init variables
while getopts ":S:m:M:T:N:O:" option; do
    case "${option}" in
        S)
            S=${OPTARG}
            ;;
        m)
            m=${OPTARG}
            ;;
        M)
			M=${OPTARG}
			;;
		T)
			T=${OPTARG}
			;;
		N)
			N=${OPTARG}
			;;
		O)
			O=${OPTARG}
			;;
        *)
            usage
            ;;
    esac
done
#testing if arguments are non-empty
shift $((OPTIND-1))
if [ -z "${S}" ] || [ -z "${m}" ] || [ -z "${M}" ] || [ -z "${N}" ] ; then
    usage
fi
#if threads empty, default value = 4
if [ -z "${T}" ]; then
	T=4
fi
#if output directory is empty, write in working directory
if [ -z "${O}" ]; then
	O="."
fi

#testing if SAM file header is correct
head -n 1 ${S} | gawk '{if($0 !~ /^@/){print "Warning : File in -S option is not a SAM format file"}}'

##Division by k-mer length
#Sample name
sample=${N}

#output directory creation
mkdir -p ${O}"bamDivision"

#Bam Division
for kmer in `seq ${m} ${M}`;
do
    #division of read depending of read length
	gawk -v l=${kmer} '{if($1 !~ /^ *@/){if(length($10)==l){print $0}}else{print $0}}' ${S} > ${O}bamDivision/${sample}.${kmer}.sam;
    #Conversion as bam
	samtools view -@ ${T} -F 268 -b ${O}bamDivision/${sample}.${kmer}.sam > ${O}bamDivision/${sample}.${kmer}.uniq.bam ;
    #sorting bam file
	samtools sort -@ ${T} ${O}bamDivision/${sample}.${kmer}.uniq.bam -o ${O}bamDivision/${sample}.${kmer}.uniq.sort.bam ;
    #index bam file
	samtools index ${O}bamDivision/${sample}.${kmer}.uniq.sort.bam;
done

#remove temporary file
rm -rf ${O}bamDivision/${sample}*.uniq.bam ${O}bamDivision/${sample}*.sam
