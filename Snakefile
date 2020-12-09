configfile: "/data/config.yaml"

# Imports
from optparse import OptionParser

#Wildcards definition
SAMPLES, = glob_wildcards("/data/fastq/{sample}.fastq.gz")
BOWTIE2 = ["1","2","3","4","rev.1","rev.2"]
HISAT2 = ["1","2","3","4","5","6","7","8"]
KMER = list(map(str,range(int(config['kmer_min']),int(config['kmer_max'])+1)))

# Change names for snamkemake workflow, depending on data given
if config['UTR'] == "yes":
    counts = "htseqcountCDS"
else:
    counts = "htseqcount"

if config['RNAseq'] == "yes":
    frag_length_S = frag_length_L = ""
else:
    frag_length_S = "." + KMER[0]
    frag_length_L = "." + KMER[0] + "-" + KMER[len(KMER)-1]


# mamba install -c bioconda -c conda-forge snakemake
# snakemake 5.26.1
rule all:
    input:
        #call of make_fastqc rule
        expand("/data/RESULTS/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        #call of quality_controls_periodicity rule
        expand("/data/RESULTS/qualitativeAnalysis/periodicity/{sample}" + frag_length_S + ".periodicity.start.CDS.-" + config['window_bf'] + "+" + config['window_af'] + ".txt", sample=SAMPLES),
        expand("/data/RESULTS/qualitativeAnalysis/periodicity/{sample}" + frag_length_S + ".periodicity.stop.CDS.-" + config['window_af'] + "+" + config['window_bf'] + ".txt", sample=SAMPLES),
        #call of quality_controls_kmerRepartition rule
        expand("/data/RESULTS/qualitativeAnalysis/kmerRepartition/{sample}.kmerRepartition.txt", sample=SAMPLES),
        #call of htseqcount_transcript_utr or htseqcount_transcript rule (depends on UTR="True"|"False" in config file)
        expand("/data/RESULTS/htseqcount_CDS/{sample}" + frag_length_L + ".no-outRNA." + counts + ".txt", sample=SAMPLES),

        # Count matrix for DESeq2
        "/data/RESULTS/Final_report.html"


# When the jobs are all done
onsuccess:
    # List of interesting logs to make the report
    logs_names = ["adapt_trimming","bowtie2_run_outRNA","run_transcriptome_hisat2","run_transcriptome_bowtie2","rpkmMoyen"]
    if config['UTR'] == "no":
        logs_names = logs_names[:-1]

    # File for the statistical report
    data_report=open("/data/RESULTS/Analysis_Report.txt","w")

    for sample in SAMPLES:
        # Data treatment report creation
        data_report.write("##################\n## NEXT SAMPLE ##\n##################\n\n" + sample + "\n")

        for log in logs_names:
            data_report.write("\n" + ("#" * (len(log)+6)) + "\n## " + log + " ##\n" + ("#" * (len(log)+6)) + "\n")
            logs_files=open("/data/logsTmp/" + sample + "_" + log + ".log","r")

            # Keep only lines of interest from cutadapt report
            i=-1
            if log=="adapt_trimming":
                if int(config['threads']) > 1:
                    lines_to_read = range(22)
                else:
                    lines_to_read = range(20)
                for position, line in enumerate(logs_files):
                    if position in lines_to_read:
                        data_report.write(line)
                    else:
                        break
            else:
                for line in logs_files:
                    data_report.write(line)

            logs_files.close()
        data_report.write("\n\n\n")
    data_report.close()

    # Removes useless directory
    shell("rm -f -r /data/RESULTS/no-outRNA/ /data/RESULTS/cutadapt/ /data/logsTmp/ ;")



# Builds the index of bowtie2 mapping for non-coding RNA
rule bowtie2_build_outRNA:
    input:
        "/data/database/" + config['fasta_outRNA']
    output:
        expand("/data/database/outRNA_bowtie2.{extb}.bt2",extb=BOWTIE2)
    params:
        outNames="/data/database/outRNA_bowtie2"
    log:
        "/data/logs/bowtie2_build_outRNA/bowtie2_build_outRNA.log"
    shell:
        # mamba install -c bioconda -c conda-forge bowtie2
        # bowtie2 2.4.2
        "bowtie2-build {input} {params.outNames} &> {log}"

# Builds the index of bowtie2 mapping for all RNA
rule bowtie2_build_transcriptome:
    input:
        "/data/database/" + config['fasta_transcriptome']
    output:
        expand("/data/database/transcriptome_bowtie2.{extb}.bt2",extb=BOWTIE2)
    params:
        outNames="/data/database/transcriptome_bowtie2"
    log:
        "/data/logs/bowtie2_build_transcriptome/bowtie2_build_transcriptome.log"
    shell:
        # mamba install -c bioconda -c conda-forge bowtie2
        # bowtie2 2.4.2
        "bowtie2-build {input} {params.outNames} &> {log}"

# Builds the index of hisat2 mapping for all RNA
rule hisat2_build_transcriptome:
    input:
        "/data/database/" + config['fasta_transcriptome']
    output:
        expand("/data/database/transcriptome_hisat2.{exth}.ht2",exth=HISAT2)
    params:
        outNames="/data/database/transcriptome_hisat2"
    log:
        "/data/logs/hisat2_build_transcriptome/hisat2_build_transcriptome.log"
    shell:
        # mamba install -c bioconda -c conda-forge hisat2
        # hisat2 2.2.1
        "hisat2-build {input} {params.outNames} &> {log}"

# Quality control of data : build of the fastqc
rule make_fastqc:
    input:
        "/data/fastq/{sample}.fastq.gz"
    output:
        "/data/RESULTS/fastqc/{sample}_fastqc.zip",
        "/data/RESULTS/fastqc/{sample}_fastqc.html"
    params:
       outdir="/data/RESULTS/fastqc/"
    log:
        "/data/logs/make_fastqc/{sample}.log"
    shell:
        # mamba install -c bioconda -c conda-forge fastqc
        # fastqc 0.11.9
        "fastqc {input} --outdir {params.outdir} 2> {log}"

# Removes/cuts potential adapters on the reads
rule adapt_trimming:
    input:
        "/data/fastq/{sample}.fastq.gz"
    output:
        "/data/RESULTS/cutadapt/{sample}.cutadapt" + frag_length_L + ".fastq.gz"
    log:
        cutadapt="/data/logs/adapt_trimming/{sample}.log",
        cutadapt_out="/data/logsTmp/{sample}_adapt_trimming.log"
    params:
        sample_names="{sample}"
    shell:
        # mamba install -c bioconda -c conda-forge cutadapt
        # cutadapt 2.10
        "cutadapt -a " + config['adapt_sequence'] + " -e 0.125 --trimmed-only --max-n=1 -m " + config['kmer_min'] + " -M " + config['kmer_max'] + " -o {output} {input} 1>> {log.cutadapt_out} 2> {log.cutadapt}"

# Mapping of non-coding RNA
rule bowtie2_run_outRNA:
    input:
        expand("/data/database/outRNA_bowtie2.{extb}.bt2",extb=BOWTIE2),
        fastq="/data/RESULTS/cutadapt/{sample}.cutadapt" + frag_length_L + ".fastq.gz"
    output:
        "/data/RESULTS/no-outRNA/{sample}" + frag_length_L + ".no-outRNA.fastq.gz"
    log:
        bt2="/data/logsTmp/{sample}_bowtie2_run_outRNA.log"
    params:
        sample_names="{sample}"
    shell:
        # mamba install -c bioconda -c conda-forge bowtie2
        # bowtie2 2.4.2
        "bowtie2 -x /data/database/outRNA_bowtie2 -p " + config['threads'] + " -U {input.fastq} --un-gz {output} > /dev/null 2>> {log.bt2} ;"
        "rm {input.fastq} ;"

# Mapping of all RNA by bowtie2 and hisat2
rule run_transcriptome:
    input:
        expand("/data/database/transcriptome_hisat2.{exth}.ht2",exth=HISAT2),
        expand("/data/database/transcriptome_bowtie2.{extb}.bt2",extb=BOWTIE2),
        fastq="/data/RESULTS/no-outRNA/{sample}" + frag_length_L + ".no-outRNA.fastq.gz"
    output:
        fastq="/data/RESULTS/no-outRNA/{sample}" + frag_length_L + ".no-outRNA.notAlign.fastq.gz",
        sam_hisat2="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".hisat2.sam",
        sam_bowtie2="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".bowtie2.sam"
    log:
        hisat2_out="/data/logsTmp/{sample}_run_transcriptome_hisat2.log",
        bowtie2_out="/data/logsTmp/{sample}_run_transcriptome_bowtie2.log"
    params:
        sample_names="{sample}"
    shell:
        # mamba install -c bioconda -c conda-forge bowtie2
        # bowtie2 2.4.2
        # mamba install -c bioconda -c conda-forge hisat2
        # hisat2 2.2.1
        "hisat2 -x /data/database/transcriptome_hisat2 -p " + config['threads'] + " -U {input.fastq} --un-gz {output.fastq} -S {output.sam_hisat2} 2>> {log.hisat2_out} ;"
        "bowtie2 -x /data/database/transcriptome_bowtie2 -p " + config['threads'] + " -U {output.fastq} -S {output.sam_bowtie2} 2>> {log.bowtie2_out};"

# Creates cram, bam and sam files
rule transcriptome_samtools:
    input:
        sam_hisat2="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".hisat2.sam",
        sam_bowtie2="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".bowtie2.sam",
        fasta="/data/database/" + config['fasta_transcriptome']
    output:
        samuniq="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".uniq.sam",
        cram="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".cram",
        crai="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".cram.crai",
        bam="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".bam",
        bai="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".bam.bai"
    log:
        egrep1="/data/logs/transcriptome_samtools/{sample}.egrep1.log",
        grep="/data/logs/transcriptome_samtools/{sample}.grep.log",
        egrep2="/data/logs/transcriptome_samtools/{sample}.egrep2.log",
        view1="/data/logs/transcriptome_samtools/{sample}.samtools_view1.log",
        view2="/data/logs/transcriptome_samtools/{sample}.samtools_view2.log",
        sort2="/data/logs/transcriptome_samtools/{sample}.samtools_sort2.log",
        view3="/data/logs/transcriptome_samtools/{sample}.samtools_view3.log",
        sort3="/data/logs/transcriptome_samtools/{sample}.samtools_sort3.log",
        index1="/data/logs/transcriptome_samtools/{sample}.samtools_index1.log",
        index2="/data/logs/transcriptome_samtools/{sample}.samtools_index2.log"
    params:
        sam="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".sam"
    shell:
        # mamba install -c bioconda -c conda-forge samtools
        # samtools 1.11
        "set +o pipefail ;"
        "egrep -i 'XM:i:2|XM:i:0|XM:i:1|@' {input.sam_hisat2} 1> {params.sam} 2> {log.egrep1} ;"
        "grep -v '@' {input.sam_bowtie2} 2> {log.grep} | egrep -i 'XM:i:2|XM:i:0|XM:i:1' 1>> {params.sam} 2> {log.egrep2} ;"
        "samtools view -@ " + config['threads'] + " -F 3332 -q 1 -h {params.sam} 1> {output.samuniq}  2> {log.view1} ;"
        "samtools view -@ " + config['threads'] + " -b {output.samuniq}  2> {log.view2} | samtools sort -@ " + config['threads'] + " -o {output.bam}  2> {log.sort2} ;"
        "samtools view -@ " + config['threads'] + " -T {input.fasta} -C {output.samuniq}  2> {log.view3} | samtools sort -@ " + config['threads'] + " -o {output.cram}  2> {log.sort3} ;"
        "samtools index {output.bam} 2> {log.index1} ;"
        "samtools index {output.cram} 2> {log.index2} ;"
        "rm {params.sam} {input.sam_hisat2} {input.sam_bowtie2};"

# Counts reads on each transcript
rule htseqcount_transcript:
    input:
        bam="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".bam",
        bai="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".bam.bai",
        gff="/data/database/" + config['gff_transcriptome'],
        CDS="/data/database/" + config['CDSlength']
    output:
        "/data/RESULTS/htseqcount_CDS/{sample}" + frag_length_L + ".no-outRNA.htseqcount.txt"
    log:
        htseqcount_CDS="/data/logs/htseqcount_transcript/{sample}.htseqcount.log"
    shell:
        # mamba install -c bioconda -c conda-forge htseq
        # htseq 0.12.4
        "set +o pipefail ;"
        "echo Genes\\\t{wildcards.sample} > {output};"
        "htseq-count -f bam -t 'CDS' -i 'ID' {input.bam} {input.gff} 2> {log.htseqcount_CDS} | head -n -5  >> {output} ;"
        "rm {input.bam} {input.bai};"

# Counts reads on each 5prime, 3prime or CDS part of each transcript
rule htseqcount_transcript_utr:
    input:
        bam="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".bam",
        bai="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".bam.bai",
        gff="/data/database/" + config['gff_transcriptome'],
        CDS="/data/database/" + config['CDSlength'],
        five_prime="/data/database/" + config['5primelength'],
        three_prime="/data/database/" + config['3primelength']
    output:
        counts_CDS="/data/RESULTS/htseqcount_CDS/{sample}" + frag_length_L + ".no-outRNA.htseqcountCDS.txt",
        counts_threeprime="/data/RESULTS/htseqcount_threeprime/{sample}" + frag_length_L + ".no-outRNA.htseqcountUTR.txt",
        counts_fiveprime="/data/RESULTS/htseqcount_fiveprime/{sample}" + frag_length_L + ".no-outRNA.htseqcountUTR.txt"
        # files_ready=touch("/data/logs/htseqcount_transcript_utr/{sample}.htseqcountUTR.mytask.done")
    log:
        view="/data/logs/htseqcount_transcript_utr/{sample}.view.log",
        htseqcount_CDS="/data/logs/htseqcount_transcript_utr/{sample}.htseqcount_CDS.log",
        join_CDS="/data/logs/htseqcount_transcript_utr/{sample}.join.log",
        awk_CDS="/data/logs/htseqcount_transcript_utr/{sample}.awk_CDS.log",
        htseqcount_3prime="/data/logs/htseqcount_transcript_utr/{sample}.htseqcount_CDS.log",
        join_3prime="/data/logs/htseqcount_transcript_utr/{sample}.join.log",
        awk_3prime="/data/logs/htseqcount_transcript_utr/{sample}.awk_CDS.log",
        htseqcount_5prime="/data/logs/htseqcount_transcript_utr/{sample}.htseqcount_CDS.log",
        join_5prime="/data/logs/htseqcount_transcript_utr/{sample}.join.log",
        awk_5prime="/data/logs/htseqcount_transcript_utr/{sample}.awk_CDS.log",
        rpkm="/data/logs/htseqcount_transcript_utr/{sample}.rpkm.log",
        rpkm_out="/data/logsTmp/{sample}_rpkmMoyen.log"
    params:
        sample_names="{sample}"
    shell:
        # mamba install -c bioconda -c conda-forge htseq
        # htseq 0.12.4
        "set +o pipefail ;"
        "totalReads=`samtools view -c {input.bam} 2> {log.view}` ;"

        "htseq-count -f bam -t 'CDS' -i 'ID' -m intersection-strict {input.bam} {input.gff} 2> {log.htseqcount_CDS} | head -n -5 > /data/RESULTS/tmpCDS.{wildcards.sample}.txt ;"
        "RPKMmoyenCDS=`join -1 1 -2 1 -o 1.1,1.2,2.2 -t $'\t' {input.CDS} /data/RESULTS/tmpCDS.{wildcards.sample}.txt 2> {log.join_CDS} | awk -F '\t' -v totalReads=$totalReads '{{print $0,$3/($2/1000*totalReads/1000000)}}' | awk 'BEGIN{{sum=0}}{{sum=sum+$4}}END{{print sum/NR}}' 2> {log.awk_CDS}` ;"
        "echo Genes\\\t{wildcards.sample} > {output.counts_CDS} ;"
        "cat /data/RESULTS/tmpCDS.{wildcards.sample}.txt >> {output.counts_CDS} ;"

        "htseq-count -f bam -t 'three_prime_UTR' -i 'ID' -m intersection-strict {input.bam} {input.gff} 2> {log.htseqcount_3prime} | head -n -5 > /data/RESULTS/tmp_threeprime.{wildcards.sample}.txt ;"
        "RPKMmoyen3prime=`join -1 1 -2 1 -o 1.1,1.2,2.2 -t $'\t' {input.three_prime} /data/RESULTS/tmp_threeprime.{wildcards.sample}.txt 2> {log.join_3prime} | awk -F '\t' -v totalReads=$totalReads '{{print $0,$3/($2/1000*totalReads/1000000)}}' | awk 'BEGIN{{sum=0}}{{sum=sum+$4}}END{{print sum/NR}}' 2> {log.awk_3prime}` ;"
        "echo Genes\\\t{wildcards.sample} > {output.counts_threeprime} ;"
        "cat /data/RESULTS/tmp_threeprime.{wildcards.sample}.txt >> {output.counts_threeprime} ;"

        "htseq-count -f bam -t 'five_prime_UTR' -i 'ID' -m intersection-strict {input.bam} {input.gff} 2> {log.htseqcount_5prime} | head -n -5 > /data/RESULTS/tmp_fiveprime.{wildcards.sample}.txt ;"
        "RPKMmoyen5prime=`join -1 1 -2 1 -o 1.1,1.2,2.2 -t $'\t' {input.five_prime} /data/RESULTS/tmp_fiveprime.{wildcards.sample}.txt 2> {log.join_5prime} | awk -F '\t' -v totalReads=$totalReads '{{print $0,$3/($2/1000*totalReads/1000000)}}' | awk 'BEGIN{{sum=0}}{{sum=sum+$4}}END{{print sum/NR}}' 2> {log.awk_5prime}` ;"
        "echo Genes\\\t{wildcards.sample} > {output.counts_fiveprime} ;"
        "cat /data/RESULTS/tmp_fiveprime.{wildcards.sample}.txt >> {output.counts_fiveprime} ;"

        "rm /data/RESULTS/tmpCDS.{wildcards.sample}.txt ;"
        "rm /data/RESULTS/tmp_threeprime.{wildcards.sample}.txt ;"
        "rm /data/RESULTS/tmp_fiveprime.{wildcards.sample}.txt ;"
        "rm {input.bam} {input.bai};"

        "echo 'RPKM moyen' 1>> {log.rpkm_out} 2> {log.rpkm} ;"
        "echo 'CDS : '$RPKMmoyenCDS 1>> {log.rpkm_out} 2>> {log.rpkm} ;"
        "echo '5prime : '$RPKMmoyen5prime 1>> {log.rpkm_out} 2>> {log.rpkm} ;"
        "echo '3prime : '$RPKMmoyen3prime 1>> {log.rpkm_out} 2>> {log.rpkm} ;"

# Quality controls :
# Divisions of SAM file according to kmer length and turns it into BAM
rule quality_controls_bamDivision:
    input:
        sam="/data/RESULTS/BAM_transcriptome/{sample}" + frag_length_L + ".uniq.sam",
        gff="/data/database/" + config['gff_transcriptome']
    output:
        bam="/data/RESULTS/qualitativeAnalysis/bamDivision/{sample}.{taille}.uniq.sort.bam",
        bai="/data/RESULTS/qualitativeAnalysis/bamDivision/{sample}.{taille}.uniq.sort.bam.bai"
    params:
        sample_names="{sample}"
    log:
        "/data/logs/quality_controls_bamDivision/{sample}.{taille}.BamDivision.log"
    shell:
        # mamba install -c bioconda -c conda-forge samtools
        # samtools 1.11
        # mamba install --name all_TRiP -c anaconda gawk
        # gawk 5.1.0
        "/TRiP/tools/BamDivision.sh -N {params.sample_names} -S {input.sam} -m " + config['kmer_min'] + " -M " + config['kmer_max'] + " -T " + config['threads'] + " -O /data/RESULTS/qualitativeAnalysis/ 2> {log};"
        "rm {input.sam};"

# Creates bed files from fasta files
rule quality_controls_bedcount:
    input:
        bam="/data/RESULTS/qualitativeAnalysis/bamDivision/{sample}.{taille}.uniq.sort.bam",
        bai="/data/RESULTS/qualitativeAnalysis/bamDivision/{sample}.{taille}.uniq.sort.bam.bai",
        fasta="/data/database/" + config['fasta_transcriptome']
    output:
        sequenceBedCount="/data/RESULTS/qualitativeAnalysis/sequenceBedCount/{sample}.{taille}.count.sequence.bed",
        kmerRepartitionBed="/data/RESULTS/qualitativeAnalysis/kmerRepartition/{sample}.{taille}.bed",
        bed="/data/RESULTS/qualitativeAnalysis/bedCount/{sample}.{taille}.count.bed",
    params:
        sample_names="{sample}"
    log:
        "/data/logs/quality_controls_bedcount/{sample}.{taille}.kmerRepartition.log"
    shell:
        # mamba install -c bioconda -c conda-forge bedtools
        # bedtools 2.29.2
        "/TRiP/tools/kmerRepartition.sh -N {params.sample_names} -F {input.fasta} -D /data/RESULTS/qualitativeAnalysis/bamDivision/ -m " + config['kmer_min'] + " -M " + config['kmer_max'] + " -O /data/RESULTS/qualitativeAnalysis/ 2> {log}"

# Outputs the number of reads on each kmer
rule quality_controls_kmerRepartition:
    input:
        "/data/RESULTS/qualitativeAnalysis/kmerRepartition/{sample}." + KMER[len(KMER)-1] + ".bed"
    output:
        "/data/RESULTS/qualitativeAnalysis/kmerRepartition/{sample}.kmerRepartition.txt"
    params:
        path="/data/RESULTS/qualitativeAnalysis/kmerRepartition/",
        sample_names="{sample}"
    log:
        wc="/data/logs/quality_controls_kmerRepartition/{sample}.wc.log",
        sed1="/data/logs/quality_controls_kmerRepartition/{sample}.sed1.log",
        awk="/data/logs/quality_controls_kmerRepartition/{sample}.awk.log",
        sed2="/data/logs/quality_controls_kmerRepartition/{sample}.sed2.log",
        head="/data/logs/quality_controls_kmerRepartition/{sample}.head.log"
    shell:
        "set +o pipefail ;"
        "wc -l {params.path}{params.sample_names}* 2> {log.wc} | sed 's/\./ /g' 2> {log.sed1} | awk -F ' ' '{{print $(NF-1),$1}}' 2> {log.awk} | sed 's/ /\t/g' 2> {log.sed2} | head -n -1 2> {log.head}  > {output} ;"

# Looks how many reads start on each base to find if there is a periodicity signal
rule quality_controls_periodicity:
    input:
        sequenceBedCount="/data/RESULTS/qualitativeAnalysis/sequenceBedCount/{sample}." + KMER[len(KMER)-1] + ".count.sequence.bed",
        kmerRepartitionBed="/data/RESULTS/qualitativeAnalysis/kmerRepartition/{sample}." + KMER[len(KMER)-1] + ".bed",
        bed="/data/RESULTS/qualitativeAnalysis/bedCount/{sample}." + KMER[len(KMER)-1] + ".count.bed",
        gff="/data/database/" + config['gff_transcriptome']
    output:
        start="/data/RESULTS/qualitativeAnalysis/periodicity/{sample}.{taille}.periodicity.start.CDS.-" + config['window_bf'] + "+" + config['window_af'] + ".txt",
        stop="/data/RESULTS/qualitativeAnalysis/periodicity/{sample}.{taille}.periodicity.stop.CDS.-" + config['window_af'] + "+" + config['window_bf'] + ".txt"
    params:
        sample_names="{sample}"
    log:
        start="/data/logs/quality_controls_periodicity/{sample}.{taille}.log",
        stop="/data/logs/quality_controls_periodicity/{sample}.{taille}.log"
    shell:
        "/TRiP/tools/periodicity.sh -N {params.sample_names} -G {input.gff} -D /data/RESULTS/qualitativeAnalysis/bedCount/ -p 'start' -t 'CDS' -m " + config['window_bf'] + " -M " + config['window_af'] + " -r 'metagene' -O /data/RESULTS/qualitativeAnalysis/ 2> {log.start} ;"
        "/TRiP/tools/periodicity.sh -N {params.sample_names} -G {input.gff} -D /data/RESULTS/qualitativeAnalysis/bedCount/ -p 'stop' -t 'CDS' -m " + config['window_af'] + " -M " + config['window_bf'] + " -r 'metagene' -O /data/RESULTS/qualitativeAnalysis/ 2> {log.stop} ;"

# Creates the row names (genes/transcript names) of the count matrix
rule count_matrix_initialization:
    input:
        ready= expand(rules.htseqcount_transcript_utr.output, sample=SAMPLES) if config['UTR']=="yes" else expand(rules.htseqcount_transcript.output, sample=SAMPLES),
        counts="/data/RESULTS/htseqcount_CDS/" + SAMPLES[0] + "" + frag_length_L + ".no-outRNA." + counts + ".txt"
    output:
        "/data/RESULTS/DESeq2/gene_list.txt"
    log:
        "/data/logs/count_matrix_initialization/cut_sort.log"
    shell:
        "set +o pipefail ;"
        "echo Genes\\\t > {output} ;"
        "cut -f1 {input.counts} | sort >> {output} 2> {log} ;"

# Adds the read counts of each sample to the matrix
rule count_matrix_creation:
    input:
        "/data/RESULTS/DESeq2/gene_list.txt"
    output:
        "/data/RESULTS/DESeq2/count_matrix.txt"
    log:
        "/data/logs/count_matrix_creation/count_matrix.log"
    params:
        tmp_file="/data/RESULTS/DESeq2/tmp.txt"
    run:
        shell("cat {input} > {output} ;")
        for sample in SAMPLES:
            shell("join --nocheck-order -t $'\t' -j 1 /data/RESULTS/htseqcount_CDS/" + sample + frag_length_L + ".no-outRNA." + counts + ".txt {output} > {params.tmp_file} 2> {log} ;")
            shell("cat {params.tmp_file} > {output} ;")
        shell("rm /data/RESULTS/DESeq2/gene_list.txt {params.tmp_file} ;")

rule DESeq2_analysis:
    input:
        counts_matrix="/data/RESULTS/DESeq2/count_matrix.txt"
    output:
        complete="/data/RESULTS/DESeq2/complete.txt",
        up="/data/RESULTS/DESeq2/up.txt",
        down="/data/RESULTS/DESeq2/down.txt",
        report="Final_report.html"
    log:
        "/data/logs/DESeq2_analysis/DESeq2_analysis.log"
    params:
        reportPath="/data/RESULTS/"
    shell:
        "Rscript -e \"rmarkdown::render('/TRiP/tools/DESeq2_analysis.Rmd', run_pandoc = FALSE, output_file={output.report}, output_dir={params.reportPath}, knit_root_dir = {params.reportPath}{output.report})\" 2> {log} ;"
        "ls /TRiP/tools/"
        "ls /data/"
        "mv /TRiP/tools/DESeq2_analysis.html /data/RESULTS/ ;"
        "ls /TRiP/"
        "ls /TRiP/tools/"
        "ls /data/"

        # "Rscript /TRiP/tools/DE_SEQ2_Analyse.R " + config['reference_condition'] + " ;"
