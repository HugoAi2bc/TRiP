# Step by step tutorial 

Welcome to the TRiP tool tutorial !  
>TRiP is designed to perform all classical steps of **ribosome profiling** (RiboSeq) data >analysis from the FastQ files to the differential expression analysis.

[TOC]

## 1) Install Docker  
First of all, Docker must be present in version 19 or higher.  
If you don’t already have it, now is the time to fix it !  
Docker Engine is available on different OS like macOS and Windows 10 through Docker Desktop and as a static binary installation for a variety of Linux platforms. All are available here : https://docs.docker.com/engine/install/  
>Tips:  
>&emsp;&emsp;&emsp;For Windows, WSL2 and Ubuntu from Microsoft store applications are needed too.  
  
## 2) Directory preparation  
TRiP does not need installation (yipee) but a precise folder architecture is required (boo).  
The first step is the project folder creation. It is named as your project and will be the volume linked to Docker.  
Then, two sub-folders and a file have to be created and completed respectively.  
> **Caution, those steps are majors for the good course of the analysis**  
> **Subfolders don’t have uppercase**    
  
Folder architecture at this step:  
Project_name  
  
### a) *fastq* subfolder  
This subfolder, as its name suggests, should contain your FastQs. These must be compressed in .gz.  
**Format of file name must be as following:**  
&emsp;&emsp;&emsp;+biological_condition_name+.<span style="color:green">replicat</span>.<span style="color:blue">fastq.gz</span>  
For example, for the first replicat of the wild-type condition, sample will be named *WT.1.fastq.gz*  
>Caution, for **Windows**, extensions can be hidden.    
   
Folder architecture at this step:  
Project_name  
└── fastq   
&emsp;&emsp;├── condA.1.fastq.gz   
&emsp;&emsp;├── condA.2.fastq.gz   
&emsp;&emsp;├── condB.1.fastq.gz   
&emsp;&emsp;└── condB.2.fastq.gz  
  
  
### b) *database* subfolder  
In this subfolder, you must put at least the following three files:  
- your genome fasta file: Whether it's the genome or the transcriptome, it must be your reference fasta file where reads will be aligned. It must, like the other files, be downloaded from the [Ensembl](https://www.ensembl.org/index.html) database.  
> Note: 
>&emsp;&emsp;For complexe genomes, transcriptome is preferable (less complexity and better analysis).   
>&emsp;&emsp;If the genome is used, the computer needs a larger RAM capacity than for transcriptome.    
   
   
- GFF3 file corresponding to the reference genome dropped.  
- out-RNA fasta: This fasta file must gather together RNA sequences you want to remove. As a rule, these are at least rRNA. You can also add mitochondrial RNA or any other RNA.  
  
If you have them, files containing each annotation length (see next paragraph) also be dropped into this folder.  
  
Folder architecture at this step:  
Project_name  
├── fastq   
│&emsp;&emsp;├── condA.1.fastq.gz   
│&emsp;&emsp;├── condA.2.fastq.gz   
│&emsp;&emsp;├── condB.1.fastq.gz   
│&emsp;&emsp;└── condB.2.fastq.gz  
└── database   
&emsp;&emsp;├── reference_transcriptome.fa  
&emsp;&emsp;├── RNA_to_remove.fa  
&emsp;&emsp;└── annotation_length.txt (if possible)  
  
### c) [config.yaml](link to configfile entre parenthèses) file  
Configfile file is used to define parameters to tell TRiP how to process your data.  
You must download it [here](lien final du config) and open it with a text editor.    
It must be carefully completed and be present in the project directory everytime you want to run TRiP.  
>Caution  
>	Spaces and quotation marks **must not be changed**! Your information must be entered in quotes     
  
#### Project name  
First and easy step, the project name ! You could (and it is recommended) use the same that your folder.  
*project_name*: principal directory name  
#### Name of database subfolder file  
You must enter the full name (**with extensions**) without the path of files added in the database subfolder previously created.   
*fasta_transcriptome*: reference transcriptome (or genome) complete name.  
*gff_transcriptome*: corresponding GFF3 name.  
*fasta_outRNA*: not interesting sequence fasta name.  
#### UTR covering option  
*UTR*: This option has to be turned on if you want to compare UTRs coverage against CDS. However, to be realized, this option requires a file with the name of each gene and the length of the associated annotation (one file by annotation: CDS, 3’-UTR and 5’-UTR).  
The full name (**with extensions**) without the path of each file has to be report in the configfile.  
*CDSlength*: complete name of file with CDS length.  
*5primelength*: complete name of file with 5’-UTR length  
*3primelength*: complete name of file with 3’-UTR length  
#### Adapters and length selection  
During the TRiP process, data is trimmed and selected depending on their length.   
*already_trimmed*: In case your data was already trimmed, you can set this option on “yes”.   
*adapt_sequence*: Else, you must specify the sequence adapters in quotes on the line below.  
You also have to define the range for read length selection. Default values are already filled.  
*kmer_min*: minimum read length  
*kmer_max*: maximum read length  
#### Statistical settings  
To be able to perform statistical analyzes, you must define a reference condition as well as your thresholds.   
*reference_condition*: it correspond to the reference biological_condition_name
We have pre-define them:  
*p-val*: defined at 0.01 to keep only genes with a high confidence.  
*logFC*: defined at 0 to keep all the genes.  
#### Window for qualitative test  
During the quality analysis, the periodicity is observed on bases around start and stop.  
>The periodicity must be calculated using a metagene profile. It provides the amount of >footprints relative to all annotated start and stop codons in a selected window.   
The window selected by default is -50/+100 nts and -100/+50 nts around start and stop codons respectively.   
*window_bf*: define your window before start and after stop  
*window_af*: define your window after star and before stop  
#### Number of threads  
Thanks to the use of Snakemake, TRiP can analyse many samples at the same time. We define that ¼ of available CPUs are necessarily requisitioned for this multiple tasks in parallel.  
As the majority of tools used in the analysis pipeline have a multithreaded option, you can choose the number of threads you allow additionally.  
*threads*: number of threads you allow  
>Caution:  
>	You could attribute maximum 4 threads.  
>	Indeed:  
>	Allow 1 thread = Remain on a ¼ of threads used  
>	Allow 2 threads = Each sample will have 2 threads at its disposal. 2/4 of the threads will therefore be used.   
>	Allow 3 threads = ¾ of CPUs are used. We advise not to go beyond so as not to saturate the computer and to be able to continue to use it.  
Folder architecture at this step:  
Project_name  
├── fastq   
│   ├── condA.1.fastq.gz   
│   ├── condA.2.fastq.gz   
│   ├── condB.1.fastq.gz   
│   └── condB.2.fastq.gz  
├── database   
│   ├── reference_transcriptome.fa  
│   ├── reference_transcriptome.gff3  
│   └── RNA_to_remove.fa  
└── configfile.yaml  

Don't forget to save file !        
  
## 3) Pull TRiP  
When the folder architecture is ready, it’s time to take a TRiP !  
First, open a terminal.  
If you never had use TRiP on your workstation, you must pull it from Docker hub.  
Copy and past the following command line:    
	`docker pull equipegst/trip`  
If you have rights problem, copy and past this command:    
	`sudo docker pull equipegst/trip`     
  
## 4) Run TRiP  
Then, or if you already have used TRiP, you can run it thanks to the following command:  
	`(sudo) docker run -v /path/to/working/directory/:/data/ equipegst/trip`  
*/path/to/working/directory/* corresponds to the project_name directory full path.   
*/data/* **must not be modified in any way**  
>Caution:  
>	All paths must start and finish with a slash “/”.   
  
   
>For Windows users  
>	Path has to start at the local disks C or D: *C:\your\path\*   
>	Path has to be composed and finished with backslashes “\”.   
>	/data/ path do not change !   
  
## 5) In case of error   
Probleme docker   
probleme de mémoire   
Et si problème autre sur un job, voir la rule et aller dans le dossier logs     
  
## 6) Understand results  
Here the project_name folder architecture after TRiP run.  
Initial folders and files are still present and highlighted in green in the tree architecture below.  
<span style="color: green">Project_name  
├── fastq/   
│   ├── condA.1.fastq.gz   
│   ├── condA.2.fastq.gz   
│   ├── condB.1.fastq.gz   
│   └── condB.2.fastq.gz  
├── database/   
│   ├── reference_transcriptome.fa</span>  
│   ├── reference_transcriptome.index.ht2  
│   ├── reference_transcriptome.index.bt2  
<span style="color: green">│   ├── reference_transcriptome.gff3  
│   ├── RNA_to_remove.fa</span>  
│   ├── RNA_to_remove.index.bt2  
<span style="color: green">│   └── annotation_length.txt (if possible)</span>  
├── logs/   
│   └── *one_log_by_step*  
├── RESULTS/  
│   ├── BAM/  
│   │   ├── *one_bam_by_sample.bam*  
│   │   └── *one_bai_by_bam.bai*  
│   ├── DESeq2/  
│   │   ├── count_matrix.txt  
│   │   ├── complete.txt  
│   │   ├── up.txt  
│   │   ├── down.txt  
│   │   └── Images/  
│   ├── fastqc/  
│   │   ├── *one_html_by_sample.html*  
│   │   └── *one_zip_by_sample.zip*  
│   ├── htseqcount_CDS/  
│   │   └── *one_file_by_sample.txt*  
│   ├── htseqcount_fiveprime/  
│   │   └── *one_file_by_sample.txt*  
│   ├── htseqcount_threeprime/  
│   │   └── *one_file_by_sample.txt*  
│   ├── qualitativeAnalysis/  
│   │   ├── bamDivision/  
│   │   │   ├── *one_bam_by_readsLength_by_sample.bam*  
│   │   │   └── *one_bai_by_bam.bai*  
│   │   ├── graphes/  
│   │   │   ├── readsLengthRepartition/  
│   │   │   │   └── *one_jpeg_by_sample.jpeg*  
│   │   │   └── periodicity/  
│   │   │       ├── *one_jpeg_by_readsLength_by_sample.start.jpeg*  
│   │   │       └── *one_jpeg_by_readsLength_by_sample.stop.jpeg*  
│   │   ├── readsLengthRepartition/  
│   │   │   └── *one_file_by_sample.txt*  
│   │   └── periodicity/  
│   │       ├── *one_file_by_readsLength_by_sample.start.txt*  
│   │       └── *one_file_by_readsLength_by_sample.stop.txt*  
│   ├── PROJECT_NAME.Analysis_report.txt  
│   ├── PROJECT_NAME.Final_report.html  
├── dag_all.svg  
<span style="color: green">└── configfile.yaml</span>  
- *logs* folder groups together all the output messages from tools used in TRiP analysis pipeline. Thus, in the event of an error, it allows you to identify the problematic step to give us feedback.   
- *RESULTS* folder contains 7 subfolders:   
	i) *BAM*: it contains a BAM file for each sample (allows visualization on tools such as IGV).  
	ii) *DESeq2*: it contains the count matrix, differential analysis tables and images related to the DESeq2 analysis.   
iii)*fastqc*: it contains raw data quality controls.   
iv) *htseqcount_CDS*:   
v) *htseqcount_fiveprime*:  
vi) *htseqcount_threeprime*:  
vii) *qualitativeAnalysis*:   
It contains also two files:  
i) *PROJECT_NAME.Analysis_report.html* gathers standard output of each analysis pipeline tool. It allows to know numbers of reads at each step a)raw reads b)reads after trimming and length selection c)after out RNA depletion d)after double alignment on the reference genome.  
ii) *PROJECT_NAME.Final_report.txt* presents all figures and explanation link with the differential analysis  
  
- *a dag file* which represents all the analysis pipeline steps with your samples.  
  
>Last big tip:  
In case that a sample is too variable against other replicats or if new samples sequencing are added to your study, you can delete/add them in the *fastq* subfolder, delete the subfolder *RESULTS/DESeq2* and the two reports presents in *RESULTS*. Run again TRiP on the same *project_name* folder and it only (re)create missing files (complete analysis for added samples, new differential analysis with all samples available in *fastq* subfolder).  
