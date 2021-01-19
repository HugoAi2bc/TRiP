################################################################################
#########################           TRiP           ################e#############
################################################################################

################################################################################
### Tool for Ribosome Profiling
### Pauline FRANCOIS ; Hugo ARBES
### December 1st, 2020
################################################################################

######
# Run this script with an interactive shell ()-i) :
# bash -i TRiP.sh
######


# Need installation de locate : $ sudo apt install locate


# conda --version;
# source /home/hugo.arbes/miniconda3/etc/profile.d/conda.sh
# MYSHELL=`basename $SHELL`;
# conda init -v $MYSHELL;

# conda activate all_TRiP;

cpu=`nproc --all`;
echo "Total CPU available : "${cpu};
cpu_use=$((${cpu}/4));
echo "Minimum CPU used for the analysis : "${cpu_use};

# conda list;
snakemake -s /TRiP/TRiP/Snakefile -j --dag -np |  dot -Tsvg > /data/dag_last-run.svg;
snakemake -s /TRiP/TRiP/Snakefile -j --dag -np --forceall |  dot -Tsvg > /data/dag_all.svg;
snakemake -s /TRiP/TRiP/Snakefile -j ${cpu_use};

# conda deactivate;


# conda activate analysis
# snakemake -s snakefile_analysis -j-np
# conda deactivate
