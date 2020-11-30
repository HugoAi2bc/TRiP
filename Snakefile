configfile: "config.yaml"
from optparse import OptionParser
#snakemake -j
#snakemake -j --dag -np --reason --forceall |  dot -Tsvg > dag.svg

# SAMPLES = ["F-rep1_S15913_all_R1_001","F-rep2_S261014_all_R1_001","M-rep1_S371115_all_R1_001","M-rep2_S481216_all_R1_001"]
# SAMPLES = ["F-rep1_S15913_all_R1_001"]

#Wildcards definition
SAMPLES, = glob_wildcards(config['workdir'] + "fastq/{sample}.fastq.gz")
BOWTIE2 = ["1","2","3","4","rev.1","rev.2"]
HISAT2 = ["1","2","3","4","5","6","7","8"]
KMER = list(map(str,range(int(config['kmer_min']),int(config['kmer_max'])+1)))

if config['UTR'] == "Yes":
    counts = "htseqcountCDS"
else:
    counts = "htseqcount"

if config['RNAseq'] == "yes":
    var_testS = var_testL = ""
else:
    var_testS = "." + KMER[0]
    var_testL = "." + KMER[0] + "-" + KMER[len(KMER)-1]


# mamba install -c bioconda -c conda-forge snakemake
# snakemake 5.26.1
rule all:
    input:
        #call of make_fastqc rule
        expand(config['workdir'] + "results/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        #call of quality_controls_periodicity rule
        expand(config['workdir'] + "results/qualitativeAnalysis/periodicity/{sample}" + var_testS + ".periodicity.start.CDS.-" + config['window_bf'] + "+" + config['window_af'] + ".txt", sample=SAMPLES),
        expand(config['workdir'] + "results/qualitativeAnalysis/periodicity/{sample}" + var_testS + ".periodicity.stop.CDS.-" + config['window_af'] + "+" + config['window_bf'] + ".txt", sample=SAMPLES),
        #call of quality_controls_kmerRepartition rule
        expand(config['workdir'] + "results/qualitativeAnalysis/kmerRepartition/{sample}.kmerRepartition.txt", sample=SAMPLES),
        #call of htseqcount_transcript_utr or htseqcount_transcript rule (depends on UTR="True"|"False" in config file)
        expand(config['workdir'] + "results/htseqcount_CDS/{sample}" + var_testL + ".no-outRNA." + counts + ".txt", sample=SAMPLES),

        # # # Files creation for SARtools
        # config['workdir'] + "results/DESeq2/target.txt",
        # config['workdir'] + "results/DESeq2/template_DESeq2.r",
        #
        # Count matrix for DESeq2
        # config['workdir'] + "results/DESeq2/count_matrix.txt"
        config['workdir'] + "results/DESeq2/complete.txt"


# When the jobs are all done
onsuccess:
    # Removes useless directory
    shell("rm -f -r " + config['workdir'] + "results/no-outRNA/ " + config['workdir'] + "results/cutadapt/;")

    # List of interesting logs to make the report
    logs_names = ["adapt_trimming","bowtie2_run_outRNA","run_transcriptome_hisat2","run_transcriptome_bowtie2","rpkmMoyen"]
    if config['UTR'] == "False":
        logs_names = logs_names[:-1]

    # File for the statistical report
    data_report=open(config['workdir'] + "resultats/Analysis_Report.txt","w")

    for sample in SAMPLES:
        # Data treatment report creation
        data_report.write("##################\n## NEXT SAMPLE ##\n##################\n\n" + sample + "\n")

        for log in logs_names:
            data_report.write("\n" + ("#" * (len(log)+6)) + "\n## " + log + " ##\n" + ("#" * (len(log)+6)) + "\n")
            logs_files=open(config['workdir'] + "logsTmp/" + sample + "_" + log + ".log","r")

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



# Builds the index of bowtie2 mapping for non-coding RNA
rule bowtie2_build_outRNA:
    input:
        config['workdir'] + "database/" + config['fasta_outRNA']
    output:
        expand(config['workdir'] + "database/outRNA_bowtie2.{extb}.bt2",extb=BOWTIE2)
    params:
        outNames=config['workdir'] + "database/outRNA_bowtie2"
    log:
        config['workdir'] + "logs/bowtie2_build_outRNA/bowtie2_build_outRNA.log"
    shell:
        # mamba install -c bioconda -c conda-forge bowtie2
        # bowtie2 2.4.2
        "bowtie2-build {input} {params.outNames} &> {log}"

# Builds the index of bowtie2 mapping for all RNA
rule bowtie2_build_transcriptome:
    input:
        config['workdir'] + "database/" + config['fasta_transcriptome']
    output:
        expand(config['workdir'] + "database/transcriptome_bowtie2.{extb}.bt2",extb=BOWTIE2)
    params:
        outNames=config['workdir'] + "database/transcriptome_bowtie2"
    log:
        config['workdir'] + "logs/bowtie2_build_transcriptome/bowtie2_build_transcriptome.log"
    shell:
        # mamba install -c bioconda -c conda-forge bowtie2
        # bowtie2 2.4.2
        "bowtie2-build {input} {params.outNames} &> {log}"

# Builds the index of hisat2 mapping for all RNA
rule hisat2_build_transcriptome:
    input:
        config['workdir'] + "database/" + config['fasta_transcriptome']
    output:
        expand(config['workdir'] + "database/transcriptome_hisat2.{exth}.ht2",exth=HISAT2)
    params:
        outNames=config['workdir'] + "database/transcriptome_hisat2"
    log:
        config['workdir'] + "logs/hisat2_build_transcriptome/hisat2_build_transcriptome.log"
    shell:
        # mamba install -c bioconda -c conda-forge hisat2
        # hisat2 2.2.1
        "hisat2-build {input} {params.outNames} &> {log}"

# Quality control of data : build of the fastqc
rule make_fastqc:
    input:
        config['workdir'] + "fastq/{sample}.fastq.gz"
    output:
        config['workdir'] + "results/fastqc/{sample}_fastqc.zip",
        config['workdir'] + "results/fastqc/{sample}_fastqc.html"
    params:
       outdir=config['workdir'] + "results/fastqc/"
    log:
        config['workdir'] + "logs/make_fastqc/{sample}.log"
    shell:
        # mamba install -c bioconda -c conda-forge fastqc
        # fastqc 0.11.9
        "fastqc {input} --outdir {params.outdir} 2> {log}"

# Removes/cuts potential adapters on the reads
rule adapt_trimming:
    input:
        config['workdir'] + "fastq/{sample}.fastq.gz"
    output:
        config['workdir'] + "results/cutadapt/{sample}.cutadapt" + var_testL + ".fastq.gz"
    log:
        cutadapt=config['workdir'] + "logs/adapt_trimming/{sample}.log",
        cutadapt_out=config['workdir'] + "logsTmp/{sample}_adapt_trimming.log"
    params:
        sample_names="{sample}"
    shell:
        # mamba install -c bioconda -c conda-forge cutadapt
        # cutadapt 2.10
        "cutadapt -a " + config['adapt_sequence'] + " -e 0.125 --trimmed-only --max-n=1 -m " + config['kmer_min'] + " -M " + config['kmer_max'] + " -o {output} {input} 1>> {log.cutadapt_out} 2> {log.cutadapt}"

# Mapping of non-coding RNA
rule bowtie2_run_outRNA:
    input:
        expand(config['workdir'] + "database/outRNA_bowtie2.{extb}.bt2",extb=BOWTIE2),
        fastq=config['workdir'] + "results/cutadapt/{sample}.cutadapt" + var_testL + ".fastq.gz"
    output:
        config['workdir'] + "results/no-outRNA/{sample}" + var_testL + ".no-outRNA.fastq.gz"
    log:
        bt2=config['workdir'] + "logsTmp/{sample}_bowtie2_run_outRNA.log"
    params:
        sample_names="{sample}"
    shell:
        # mamba install -c bioconda -c conda-forge bowtie2
        # bowtie2 2.4.2
        "bowtie2 -x " + config['workdir'] + "database/outRNA_bowtie2 -p " + config['threads'] + " -U {input.fastq} --un-gz {output} > /dev/null 2>> {log.bt2} ;"
        "rm {input.fastq} ;"

# Mapping of all RNA by bowtie2 and hisat2
rule run_transcriptome:
    input:
        expand(config['workdir'] + "database/transcriptome_hisat2.{exth}.ht2",exth=HISAT2),
        expand(config['workdir'] + "database/transcriptome_bowtie2.{extb}.bt2",extb=BOWTIE2),
        fastq=config['workdir'] + "results/no-outRNA/{sample}" + var_testL + ".no-outRNA.fastq.gz"
    output:
        fastq=config['workdir'] + "results/no-outRNA/{sample}" + var_testL + ".no-outRNA.notAlign.fastq.gz",
        sam_hisat2=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".hisat2.sam",
        sam_bowtie2=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".bowtie2.sam"
    log:
        hisat2_out=config['workdir'] + "logsTmp/{sample}_run_transcriptome_hisat2.log",
        bowtie2_out=config['workdir'] + "logsTmp/{sample}_run_transcriptome_bowtie2.log"
    params:
        sample_names="{sample}"
    shell:
        # mamba install -c bioconda -c conda-forge bowtie2
        # bowtie2 2.4.2
        # mamba install -c bioconda -c conda-forge hisat2
        # hisat2 2.2.1
        "hisat2 -x " + config['workdir'] + "database/transcriptome_hisat2 -p " + config['threads'] + " -U {input.fastq} --un-gz {output.fastq} -S {output.sam_hisat2} 2>> {log.hisat2_out} ;"
        "bowtie2 -x " + config['workdir'] + "database/transcriptome_bowtie2 -p " + config['threads'] + " -U {output.fastq} -S {output.sam_bowtie2} 2>> {log.bowtie2_out};"

# Creates cram, bam and sam files
rule transcriptome_samtools:
    input:
        sam_hisat2=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".hisat2.sam",
        sam_bowtie2=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".bowtie2.sam",
        fasta=config['workdir'] + "database/" + config['fasta_transcriptome']
    output:
        samuniq=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".uniq.sam",
        cram=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".cram",
        crai=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".cram.crai",
        bam=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".bam",
        bai=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".bam.bai"
    log:
        egrep1=config['workdir'] + "logs/transcriptome_samtools/{sample}.egrep1.log",
        grep=config['workdir'] + "logs/transcriptome_samtools/{sample}.grep.log",
        egrep2=config['workdir'] + "logs/transcriptome_samtools/{sample}.egrep2.log",
        view1=config['workdir'] + "logs/transcriptome_samtools/{sample}.samtools_view1.log",
        view2=config['workdir'] + "logs/transcriptome_samtools/{sample}.samtools_view2.log",
        sort2=config['workdir'] + "logs/transcriptome_samtools/{sample}.samtools_sort2.log",
        view3=config['workdir'] + "logs/transcriptome_samtools/{sample}.samtools_view3.log",
        sort3=config['workdir'] + "logs/transcriptome_samtools/{sample}.samtools_sort3.log",
        index1=config['workdir'] + "logs/transcriptome_samtools/{sample}.samtools_index1.log",
        index2=config['workdir'] + "logs/transcriptome_samtools/{sample}.samtools_index2.log"
    params:
        sam=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".sam"
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
        bam=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".bam",
        bai=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".bam.bai",
        gff=config['workdir'] + "database/" + config['gff_transcriptome'],
        CDS=config['workdir'] + "database/" + config['CDSlength']
    output:
        config['workdir'] + "results/htseqcount_CDS/{sample}" + var_testL + ".no-outRNA.htseqcount.txt"
        # files_ready=touch(config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.htseqcount.mytask.done")
    log:
        htseqcount_CDS=config['workdir'] + "logs/htseqcount_transcript/{sample}.htseqcount.log"
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
        bam=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".bam",
        bai=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".bam.bai",
        gff=config['workdir'] + "database/" + config['gff_transcriptome'],
        CDS=config['workdir'] + "database/" + config['CDSlength'],
        five_prime=config['workdir'] + "database/" + config['5primelength'],
        three_prime=config['workdir'] + "database/" + config['3primelength']
    output:
        counts_CDS=config['workdir'] + "results/htseqcount_CDS/{sample}" + var_testL + ".no-outRNA.htseqcountCDS.txt",
        counts_threeprime=config['workdir'] + "results/htseqcount_threeprime/{sample}" + var_testL + ".no-outRNA.htseqcountUTR.txt",
        counts_fiveprime=config['workdir'] + "results/htseqcount_fiveprime/{sample}" + var_testL + ".no-outRNA.htseqcountUTR.txt"
        # files_ready=touch(config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.htseqcountUTR.mytask.done")
    log:
        view=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.view.log",
        htseqcount_CDS=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.htseqcount_CDS.log",
        join_CDS=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.join.log",
        awk_CDS=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.awk_CDS.log",
        htseqcount_3prime=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.htseqcount_CDS.log",
        join_3prime=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.join.log",
        awk_3prime=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.awk_CDS.log",
        htseqcount_5prime=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.htseqcount_CDS.log",
        join_5prime=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.join.log",
        awk_5prime=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.awk_CDS.log",
        rpkm=config['workdir'] + "logs/htseqcount_transcript_utr/{sample}.rpkm.log",
        rpkm_out=config['workdir'] + "logsTmp/{sample}_rpkmMoyen.log"
    params:
        sample_names="{sample}"
    shell:
        # mamba install -c bioconda -c conda-forge htseq
        # htseq 0.12.4
        "set +o pipefail ;"
        "totalReads=`samtools view -c {input.bam} 2> {log.view}` ;"

        "htseq-count -f bam -t 'CDS' -i 'ID' -m intersection-strict {input.bam} {input.gff} 2> {log.htseqcount_CDS} | head -n -5 > " + config['workdir'] + "results/tmpCDS.{wildcards.sample}.txt ;"
        "RPKMmoyenCDS=`join -1 1 -2 1 -o 1.1,1.2,2.2 -t $'\t' {input.CDS} results/tmpCDS.{wildcards.sample}.txt 2> {log.join_CDS} | awk -F '\t' -v totalReads=$totalReads '{{print $0,$3/($2/1000*totalReads/1000000)}}' | awk 'BEGIN{{sum=0}}{{sum=sum+$4}}END{{print sum/NR}}' 2> {log.awk_CDS}` ;"
        "echo Genes\\\t{wildcards.sample} > {output.counts_CDS} ;"
        "cat " + config['workdir'] + "results/tmpCDS.{wildcards.sample}.txt >> {output.counts_CDS} ;"

        "htseq-count -f bam -t 'three_prime_UTR' -i 'ID' -m intersection-strict {input.bam} {input.gff} 2> {log.htseqcount_3prime} | head -n -5 > " + config['workdir'] + "results/tmp_threeprime.{wildcards.sample}.txt ;"
        "RPKMmoyen3prime=`join -1 1 -2 1 -o 1.1,1.2,2.2 -t $'\t' {input.three_prime} results/tmp_threeprime.{wildcards.sample}.txt 2> {log.join_3prime} | awk -F '\t' -v totalReads=$totalReads '{{print $0,$3/($2/1000*totalReads/1000000)}}' | awk 'BEGIN{{sum=0}}{{sum=sum+$4}}END{{print sum/NR}}' 2> {log.awk_3prime}` ;"
        "echo Genes\\\t{wildcards.sample} > {output.counts_threeprime} ;"
        "cat " + config['workdir'] + "results/tmp_threeprime.{wildcards.sample}.txt >> {output.counts_threeprime} ;"

        "htseq-count -f bam -t 'five_prime_UTR' -i 'ID' -m intersection-strict {input.bam} {input.gff} 2> {log.htseqcount_5prime} | head -n -5 > " + config['workdir'] + "results/tmp_fiveprime.{wildcards.sample}.txt ;"
        "RPKMmoyen5prime=`join -1 1 -2 1 -o 1.1,1.2,2.2 -t $'\t' {input.five_prime} results/tmp_fiveprime.{wildcards.sample}.txt 2> {log.join_5prime} | awk -F '\t' -v totalReads=$totalReads '{{print $0,$3/($2/1000*totalReads/1000000)}}' | awk 'BEGIN{{sum=0}}{{sum=sum+$4}}END{{print sum/NR}}' 2> {log.awk_5prime}` ;"
        "echo Genes\\\t{wildcards.sample} > {output.counts_fiveprime} ;"
        "cat " + config['workdir'] + "results/tmp_fiveprime.{wildcards.sample}.txt >> {output.counts_fiveprime} ;"

        "rm " + config['workdir'] + "results/tmpCDS.{wildcards.sample}.txt ;"
        "rm " + config['workdir'] + "results/tmp_threeprime.{wildcards.sample}.txt ;"
        "rm " + config['workdir'] + "results/tmp_fiveprime.{wildcards.sample}.txt ;"
        "rm {input.bam} {input.bai};"

        "echo 'RPKM moyen' 1>> {log.rpkm_out} 2> {log.rpkm} ;"
        "echo 'CDS : '$RPKMmoyenCDS 1>> {log.rpkm_out} 2>> {log.rpkm} ;"
        "echo '5prime : '$RPKMmoyen5prime 1>> {log.rpkm_out} 2>> {log.rpkm} ;"
        "echo '3prime : '$RPKMmoyen3prime 1>> {log.rpkm_out} 2>> {log.rpkm} ;"

# Quality control :
# Divisions of SAM file according to kmer length and turns it into BAM
rule quality_controls_bamDivision:
    input:
        sam=config['workdir'] + "results/BAM_transcriptome/{sample}" + var_testL + ".uniq.sam",
        # fasta=config['workdir'] + "database/" + config['fasta_transcriptome'],
        gff=config['workdir'] + "database/" + config['gff_transcriptome']
    output:
        bam=config['workdir'] + "results/qualitativeAnalysis/bamDivision/{sample}.{taille}.uniq.sort.bam",
        bai=config['workdir'] + "results/qualitativeAnalysis/bamDivision/{sample}.{taille}.uniq.sort.bam.bai"
    params:
        sample_names="{sample}"
    log:
        config['workdir'] + "logs/quality_controls_bamDivision/{sample}.{taille}.BamDivision.log"
    shell:
        # mamba install -c bioconda -c conda-forge samtools
        # samtools 1.11
        "./tools/BamDivision.sh -N {params.sample_names} -S {input.sam} -m " + config['kmer_min'] + " -M " + config['kmer_max'] + " -T " + config['threads'] + " -O " + config['workdir'] + "results/qualitativeAnalysis/ 2> {log};"
        # "/tools/BamDivision.sh -N {params.sample_names} -S {input.sam} -m " + config['kmer_min'] + " -M " + config['kmer_max'] + " -T " + config['threads'] + " -O " + config['workdir'] + "results/qualitativeAnalysis/ 2> {log}"
        "rm {input.sam};"

# Creates bed files from fasta files
rule quality_controls_bedcount:
    input:
        bam=config['workdir'] + "results/qualitativeAnalysis/bamDivision/{sample}.{taille}.uniq.sort.bam",
        bai=config['workdir'] + "results/qualitativeAnalysis/bamDivision/{sample}.{taille}.uniq.sort.bam.bai",
        fasta=config['workdir'] + "database/" + config['fasta_transcriptome']
    output:
        sequenceBedCount=config['workdir'] + "results/qualitativeAnalysis/sequenceBedCount/{sample}.{taille}.count.sequence.bed",
        kmerRepartitionBed=config['workdir'] + "results/qualitativeAnalysis/kmerRepartition/{sample}.{taille}.bed",
        bed=config['workdir'] + "results/qualitativeAnalysis/bedCount/{sample}.{taille}.count.bed",
    params:
        sample_names="{sample}"
    log:
        config['workdir'] + "logs/quality_controls_bedcount/{sample}.{taille}.kmerRepartition.log"
    shell:
        # mamba install -c bioconda -c conda-forge bedtools
        # bedtools 2.29.2
        "./tools/kmerRepartition.sh -N {params.sample_names} -F {input.fasta} -D " + config['workdir'] + "results/qualitativeAnalysis/bamDivision/ -m " + config['kmer_min'] + " -M " + config['kmer_max'] + " -O " + config['workdir'] + "results/qualitativeAnalysis/ 2> {log}"
        # "/tools/kmerRepartition.sh -N {params.sample_names} -F {input.fasta} -D " + config['workdir'] + "results/qualitativeAnalysis/bamDivision/ -m " + config['kmer_min'] + " -M " + config['kmer_max'] + " -O " + config['workdir'] + "results/qualitativeAnalysis/ 2> {log}"

# Outputs the number of reads on each kmer
rule quality_controls_kmerRepartition:
    input:
        config['workdir'] + "results/qualitativeAnalysis/kmerRepartition/{sample}." + KMER[len(KMER)-1] + ".bed"
    output:
        config['workdir'] + "results/qualitativeAnalysis/kmerRepartition/{sample}.kmerRepartition.txt"
    params:
        path=config['workdir'] + "results/qualitativeAnalysis/kmerRepartition/",
        sample_names="{sample}"
    log:
        wc=config['workdir'] + "logs/quality_controls_kmerRepartition/{sample}.wc.log",
        sed1=config['workdir'] + "logs/quality_controls_kmerRepartition/{sample}.sed1.log",
        awk=config['workdir'] + "logs/quality_controls_kmerRepartition/{sample}.awk.log",
        sed2=config['workdir'] + "logs/quality_controls_kmerRepartition/{sample}.sed2.log",
        head=config['workdir'] + "logs/quality_controls_kmerRepartition/{sample}.head.log"
    shell:
        "set +o pipefail ;"
        "wc -l {params.path}{params.sample_names}* 2> {log.wc} | sed 's/\./ /g' 2> {log.sed1} | awk -F ' ' '{{print $(NF-1),$1}}' 2> {log.awk} | sed 's/ /\t/g' 2> {log.sed2} | head -n -1 2> {log.head}  > {output} ;"

# Looks how many reads start on each base to find if there is a periodicity signal
rule quality_controls_periodicity:
    input:
        sequenceBedCount=config['workdir'] + "results/qualitativeAnalysis/sequenceBedCount/{sample}." + KMER[len(KMER)-1] + ".count.sequence.bed",
        kmerRepartitionBed=config['workdir'] + "results/qualitativeAnalysis/kmerRepartition/{sample}." + KMER[len(KMER)-1] + ".bed",
        bed=config['workdir'] + "results/qualitativeAnalysis/bedCount/{sample}." + KMER[len(KMER)-1] + ".count.bed",
        gff=config['workdir'] + "database/" + config['gff_transcriptome']
    output:
        start=config['workdir'] + "results/qualitativeAnalysis/periodicity/{sample}.{taille}.periodicity.start.CDS.-" + config['window_bf'] + "+" + config['window_af'] + ".txt",
        stop=config['workdir'] + "results/qualitativeAnalysis/periodicity/{sample}.{taille}.periodicity.stop.CDS.-" + config['window_af'] + "+" + config['window_bf'] + ".txt"
    params:
        sample_names="{sample}"
    log:
        start=config['workdir'] + "logs/quality_controls_periodicity/{sample}.{taille}.log",
        stop=config['workdir'] + "logs/quality_controls_periodicity/{sample}.{taille}.log"
    shell:
        "./tools/periodicity.sh -N {params.sample_names} -G {input.gff} -D " + config['workdir'] + "results/qualitativeAnalysis/bedCount/ -p 'start' -t 'CDS' -m " + config['window_bf'] + " -M " + config['window_af'] + " -r 'metagene' -O " + config['workdir'] + "results/qualitativeAnalysis/ 2> {log.start} ;"
        # "/tools/periodicity.sh -N {params.sample_names} -G {input.gff} -D " + config['workdir'] + "results/qualitativeAnalysis/bedCount/ -p 'start' -t 'CDS' -m " + config['window_bf'] + " -M " + config['window_af'] + " -r 'metagene' -O " + config['workdir'] + "results/qualitativeAnalysis/ 2> {log.start} ;"
        "./tools/periodicity.sh -N {params.sample_names} -G {input.gff} -D " + config['workdir'] + "results/qualitativeAnalysis/bedCount/ -p 'stop' -t 'CDS' -m " + config['window_af'] + " -M " + config['window_bf'] + " -r 'metagene' -O " + config['workdir'] + "results/qualitativeAnalysis/ 2> {log.stop} ;"
        # "/tools/periodicity.sh -N {params.sample_names} -G {input.gff} -D " + config['workdir'] + "results/qualitativeAnalysis/bedCount/ -p 'stop' -t 'CDS' -m " + config['window_af'] + " -M " + config['window_bf'] + " -r 'metagene' -O " + config['workdir'] + "results/qualitativeAnalysis/ 2> {log.stop}"

# Creates the row names (genes/transcript names) of the count matrix
rule count_matrix_initialization:
    input:
        ready= expand(rules.htseqcount_transcript_utr.output, sample=SAMPLES) if config['UTR']=="True" else expand(rules.htseqcount_transcript.output, sample=SAMPLES),
        # ready=config['workdir'] + "logs/htseqcount_transcript_utr/" + SAMPLES[0] + "." + counts + ".mytask.done",
        counts=config['workdir'] + "results/htseqcount_CDS/" + SAMPLES[0] + "" + var_testL + ".no-outRNA." + counts + ".txt"
    output:
        config['workdir'] + "results/DESeq2/gene_list.txt"
    log:
        config['workdir'] + "logs/count_matrix_initialization/cut_sort.log"
    shell:
        "set +o pipefail ;"
        "echo Genes\\\t > {output} ;"
        "cut -f1 {input.counts} | sort >> {output} 2> {log} ;"
        # "rm {input.ready} ;"

# Adds the read counts of each sample to the matrix
rule count_matrix_creation:
    input:
        # rules.count_matrix_initialization.output
        config['workdir'] + "results/DESeq2/gene_list.txt"
    output:
        config['workdir'] + "results/DESeq2/count_matrix.txt"
    log:
        config['workdir'] + "logs/count_matrix_creation/count_matrix.log"
    params:
        # sample=config['workdir'] + "results/htseqcount_CDS/{sample}" + var_testL + ".no-outRNA." + counts + ".txt"
        tmp_file=config['workdir'] + "results/DESeq2/tmp.txt"
    run:
        shell("cat {input} > {output} ;")
        for sample in SAMPLES:
            shell("join --nocheck-order -t $'\t' -j 1 " + config['workdir'] + "results/htseqcount_CDS/" + sample + "" + var_testL + ".no-outRNA." + counts + ".txt {output} > {params.tmp_file} 2> {log} ;")
            shell("cat {params.tmp_file} > {output} ;")
        shell("rm " + config['workdir'] + "results/DESeq2/gene_list.txt {params.tmp_file} ;")

rule DESeq2_target_creation:
    output:
        # "/home/hugo.arbes/Hugo/GST/snakemake_pipeline/results/DESeq2/target.txt"
        config['workdir'] + "results/DESeq2/target.txt"
    log:
        config['workdir'] + "logs/DESeq2_target_creation/DESeq2_target_creation.log"
    shell:
        "python " + config['workdir'] + "target_creation.py -n \"" + ' '.join(SAMPLES) + "\" -r " + config['reference_condition'] + " -m " + KMER[0] + " -M " + KMER[len(KMER)-1] + " -p " + config['workdir'] + "results/DESeq2/ 2> {log} ;"
        #"python /tools/template_creation.py -n \"" + ' '.join(SAMPLES) + "\" -r " + config['reference_condition'] + " -m " + KMER[0] + " -M " + KMER[len(KMER)-1]

rule DESeq2_template_creation:
    output:
        # "/home/hugo.arbes/Hugo/GST/snakemake_pipeline/results/DESeq2/template_DESeq2.r"
        config['workdir'] + "results/DESeq2/template_DESeq2.r"
    log:
        config['workdir'] + "logs/DESeq2_template_creation/DESeq2_template_creation.log"
    shell:
        "python " + config['workdir'] + "template_creation.py -r " + config['reference_condition'] + " -n " + config['project_name'] + " -p " + config['workdir'] + " 2> {log} ;"
        # "python /tools/template_creation.py -r " + config['reference_condition'] + " -n " + config['project_name']



rule DESeq2_analysis:
    input:
        counts_matrix=config['workdir'] + "results/DESeq2/count_matrix.txt"
    output:
        config['workdir'] + "results/DESeq2/complete.txt",
        config['workdir'] + "results/DESeq2/up.txt",
        config['workdir'] + "results/DESeq2/down.txt"
    shell:
        "Rscript ./tools/DE_SEQ2_Analyse.R " + config['reference_condition'] + " ;"
        "mv Rplots.pdf ./results/DESeq2/ ;"



# rule DESeq2_running:
#     input:
#         # template=rules.DESeq2_template_creation.output,
#         # target=rules.DESeq2_target_creation.output,
#         # count=rules.htseqcount_transcript.output
#         target=config['workdir'] + "results/DESeq2/target.txt",
#         template=config['workdir'] + "results/DESeq2/template_DESeq2.r",
#         count=config['workdir'] + "results/htseq-count_CDS/" + SAMPLES[len(SAMPLES)-1] + "" + var_testL + ".no-outRNA." + counts + ".txt"
#     output:
#         config['project_name'] + "_report.html",
#         config['project_name'] + ".Rdata"
#     singularity:
#         "docker://genomicpariscentre/sartools"
#     run:
#         R("source({input.template})")




# rule TOUT_FOUTRE_EN_HTML:
#     input:
#         config['workdir'] + "Analysis_Report.txt"
#     output:
#         config['workdir'] + "Analysis_Report.html"
#     script:
#         script_qui_dechire_pour_faire_du_html.sh
#     # shell:
#         # "script_qui_dechire_pour_faire_du_html.sh"



# rule DESeq2_editing:
#     input:
#         CDS_count=expand(config['workdir'] + "results/htseq-count_CDS/{sample}" + var_testL + ".no-outRNA." + counts + ".txt", sample=SAMPLES) if config['UTR']=="True" else expand(config['workdir'] + "results/htseq-count_CDS/{sample}" + var_testL + ".no-outRNA.htseqcount.txt", sample=SAMPLES)
#     output:
#         target=config['workdir'] + "results/DESeq2/target.txt",
#         template=config['workdir'] + "results/DESeq2/template_script_DESeq2.r"
#     shell:
#         "echo 'label    file    group   replicat' > {output.target} ;"




# rule DESeq2_running:
#
