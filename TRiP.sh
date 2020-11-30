#! /bin/bash

################################################################################
#########################           TRiP           ################e#############
################################################################################

################################################################################
### Tool for Ribosome Profiling
### Pauline FRANCOIS ; Hugo ARBES
### October 14th, 2020
################################################################################

######
# Run this script with an interactive shell ()-i) :
# bash -i TRiP.sh
######


# Need installation de locate : $ sudo apt install locate


conda --version;
# source /home/hugo.arbes/miniconda3/etc/profile.d/conda.sh
# MYSHELL=`basename $SHELL`;
# conda init -v $MYSHELL;

# conda activate all_TRiP;

# conda list;
snakemake -s TRiP/Snakefile -j --dag -np |  dot -Tsvg > dag_last-run.svg;
snakemake -s TRiP/Snakefile -j --dag -np --forceall |  dot -Tsvg > dag_all.svg;
snakemake -s TRiP/Snakefile -j -np;

# conda deactivate;


# conda activate analysis
# snakemake -s snakefile_analysis -j-np
# conda deactivate
