BEGIN{print "sampleName\treadLengths\ttotalReads\twrittenReads\twrittenReadsPercentage\tadaptatorSequence\treadsAfterDepletion\treadsRemovePercentage\treadsAligns\tRPKMcds\tRPKM5prime\tRPKM3prime"}
{
    if($2=="NEXT"){getline;getline;getline;printf $1"\t"}
    if($1=="Command" && $2=="line"){printf $11"-"$13"\t"}
    if($1=="Total" && $2=="reads"){printf $4"\t"}
    if($1=="Reads" && $2=="written"){printf $5"\t"$6"\t"}
    if($1=="Sequence:"){printf $2"\t"}
    if($2=="bowtie2_run_outRNA"){getline;getline;getline;getline;printf $1"\t";getline;printf $2"\t"}
    if($2=="run_transcriptome_hisat2"){getline;getline;getline;getline;getline;readTranscript=$1"\t"}
    if($2=="run_transcriptome_bowtie2"){getline;getline;getline;getline;getline;readTranscript+=$1;printf readTranscript"\t"}
    if($2=="rpkmMoyen"){getline;getline;getline;printf $3"\t";getline;printf $3"\t";getline;printf $3"\n"}
}
END{}

#sampleName = key

#Sous dico trimming
#adaptateur
# tailles conservées
#total reads processed
#reads written

#sous dico depletion
#pourcentage de contamination
#nombre de reads restant

#sous dico alignement
#nombre de reads alignés par hisat et bowtie2 sur le transcriptome

#sous dico rpkm
#rpkm moyen CDS 5prime 3prime
