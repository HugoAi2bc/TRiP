#! /bin/bash

#####
# Save genes, padj et logFC after DESeq2 analysis
#####


awk '{print $1"\t"$13"\t"$14"\t"$16"\t"$17;}' ./results/DESeq2/up.txt > ./results/DESeq2/tmp.txt
tail -n +2 ./results/DESeq2/tmp.txt > ./results/DESeq2/test_gene_logFC_pval_padj_UP.txt

sed -r -e 's/:.+:CDS\t/\t/g' ./results/DESeq2/test_gene_logFC_pval_padj_UP.txt > ./results/DESeq2/test_gene_logFC_pval_padj_ENSGonly_UP.txt

head -1 ./results/DESeq2/test_gene_logFC_pval_padj_UP.txt > ./results/DESeq2/test_gene_logFC_pval_padj_gene_names_UP.txt
# FAUT TRIER POUR LE JOIN !!!
join -2 2 -o 2.3,1.2,1.3,1.4,1.5 --nocheck-order ./results/DESeq2/test_gene_logFC_pval_padj_ENSGonly_UP.txt ./correspondancetranscript_GeneID_GeneName.sort.txt >>\
./results/DESeq2/test_gene_logFC_pval_padj_gene_names_UP.txt

sed -i 's/ /\t/g' ./results/DESeq2/test_gene_logFC_pval_padj_gene_names_UP.txt

tail -n +2 ./results/DESeq2/test_gene_logFC_pval_padj_gene_names_UP.txt | cut -f1 > ./results/DESeq2/test_gene_names_UP.txt



# Take infos of : Gene logFC lfcSE pval padj
awk '{print $1"\t"$13"\t"$14"\t"$16"\t"$17;}' ./results/DESeq2/up.txt > ./results/DESeq2/tmp.txt
# Sort on padj
(head -n 1 ./results/DESeq2/tmp.txt && tail -n +2 ./results/DESeq2/tmp.txt | sort -k5) > ./results/DESeq2/gene_logFC_pval_padj_UP.txt
rm ./results/DESeq2/tmp.txt
awk '{print $1"\t"$13"\t"$14"\t"$16"\t"$17;}' ./results/DESeq2/down.txt > ./results/DESeq2/tmp.txt
(head -n 1 ./results/DESeq2/tmp.txt && tail -n +2 ./results/DESeq2/tmp.txt | sort -k5) > ./results/DESeq2/gene_logFC_pval_padj_DOWN.txt
rm ./results/DESeq2/tmp.txt
awk '{print $1"\t"$13"\t"$14"\t"$16"\t"$17;}' ./results/DESeq2/complete.txt > ./results/DESeq2/tmp.txt
(head -n 1 ./results/DESeq2/tmp.txt && tail -n +2 ./results/DESeq2/tmp.txt | sort -k5) > ./results/DESeq2/gene_logFC_pval_padj_ALL.txt
rm ./results/DESeq2/tmp.txt

####
# find a gene in first table
# grep GENE_ID /media/partage/bioinfo/Shwachman-Diamond/run/results/DESeq2/FminusvsMminus.complete.txt
# head -1 /media/partage/bioinfo/Shwachman-Diamond/run/results/DESeq2/FminusvsMminus.complete.txt | cut -f 1,13,14,15,16; grep ENSG00000153560.11 /media/partage/bioinfo/Shwachman-Diamond/run/results/DESeq2/FminusvsMminus.complete.txt | cut -f 1,13,14,15,16
####

# awk '{print $1;}' ./results/DESeq2/gene_list_FC-pval_UP.txt | tail -n +2 > ./results/DESeq2/gene_list_UP.txt
# awk '{print $1;}' ./results/DESeq2/gene_list_FC-pval_DOWN.txt | tail -n +2 > ./results/DESeq2/gene_list_DOWN.txt
# awk '{print $1;}' ./results/DESeq2/gene_list_FC-pval_ALL.txt | tail -n +2 > ./results/DESeq2/gene_list_ALL.txt
#
# cut -f1 -d ":" ./results/DESeq2/gene_list_UP.txt | sort > ./results/DESeq2/gene_list_gene-IDs_UP.txt
# cut -f1 -d ":" ./results/DESeq2/gene_list_DOWN.txt | sort > ./results/DESeq2/gene_list_gene-IDs_DOWN.txt
# cut -f1 -d ":" ./results/DESeq2/gene_list_ALL.txt | sort > ./results/DESeq2/gene_list_gene-IDs_ALL.txt
#
# rm ./results/DESeq2/gene_list_UP.txt ./results/DESeq2/gene_list_DOWN.txt ./results/DESeq2/gene_list_ALL.txt
#
# # Selection des genes d'interet en fonction de la pval a la main, puis:
# join -2 2 -o 1.1,2.3 ./results/DESeq2/gene_list_gene-IDs_UP.txt ./correspondancetranscript_GeneID_GeneName.sort.txt |\
# cut -f2 -d " " > ./results/DESeq2/gene_list_gene-names_pval0.001_UP.txt
# join -2 2 -o 1.1,2.3 ./results/DESeq2/gene_list_gene-IDs_DOWN.txt ./correspondancetranscript_GeneID_GeneName.sort.txt |\
# cut -f2 -d " " > ./results/DESeq2/gene_list_gene-names_pval0.001_DOWN.txt
# join -2 2 -o 1.1,2.3 ./results/DESeq2/gene_list_gene-IDs_ALL.txt ./correspondancetranscript_GeneID_GeneName.sort.txt |\
# cut -f2 -d " " > ./results/DESeq2/gene_list_gene-names_pval0.001_ALL.txt
#
# rm ./results/DESeq2/gene_list_gene-IDs_UP.txt ./results/DESeq2/gene_list_gene-IDs_DOWN.txt ./results/DESeq2/gene_list_gene-IDs_ALL.txt
# rm ./results/DESeq2/gene_list_gene-names_pval0.001_UP.txt ./results/DESeq2/gene_list_gene-names_pval0.001_DOWN.txt ./results/DESeq2/gene_list_gene-names_pval0.001_ALL.txt


tail -n +2 ./results/DESeq2/gene_logFC_pval_padj_UP.txt | sed -r -e 's/:.+:CDS\t/\t/g' | sort > ./results/DESeq2/gene_logFC_pval_padj_ENSGonly_UP.txt
tail -n +2 ./results/DESeq2/gene_logFC_pval_padj_DOWN.txt | sed -r -e 's/:.+:CDS\t/\t/g' | sort > ./results/DESeq2/gene_logFC_pval_padj_ENSGonly_DOWN.txt
tail -n +2 ./results/DESeq2/gene_logFC_pval_padj_ALL.txt | sed -r -e 's/:.+:CDS\t/\t/g' | sort > ./results/DESeq2/gene_logFC_pval_padj_ENSGonly_ALL.txt

head -1 ./results/DESeq2/gene_logFC_pval_padj_UP.txt > ./results/DESeq2/gene_logFC_pval_padj_gene_names_UP.txt
join -2 2 -o 2.3,1.2,1.3,1.4,1.5 ./results/DESeq2/gene_logFC_pval_padj_ENSGonly_UP.txt ./correspondancetranscript_GeneID_GeneName.sort.txt >>\
./results/DESeq2/gene_logFC_pval_padj_gene_names_UP.txt
head -1 ./results/DESeq2/gene_logFC_pval_padj_DOWN.txt > ./results/DESeq2/gene_logFC_pval_padj_gene_names_DOWN.txt
join -2 2 -o 2.3,1.2,1.3,1.4,1.5 ./results/DESeq2/gene_logFC_pval_padj_ENSGonly_DOWN.txt ./correspondancetranscript_GeneID_GeneName.sort.txt >>\
./results/DESeq2/gene_logFC_pval_padj_gene_names_DOWN.txt
head -1 ./results/DESeq2/gene_logFC_pval_padj_ALL.txt > ./results/DESeq2/gene_logFC_pval_padj_gene_names_ALL.txt
join -2 2 -o 2.3,1.2,1.3,1.4,1.5 ./results/DESeq2/gene_logFC_pval_padj_ENSGonly_ALL.txt ./correspondancetranscript_GeneID_GeneName.sort.txt >>\
./results/DESeq2/gene_logFC_pval_padj_gene_names_ALL.txt
rm ./results/DESeq2/gene_logFC_pval_padj_UP.txt ./results/DESeq2/gene_logFC_pval_padj_DOWN.txt ./results/DESeq2/gene_logFC_pval_padj_ALL.txt

sed -i 's/ /\t/g' ./results/DESeq2/gene_logFC_pval_padj_gene_names_UP.txt
sed -i 's/ /\t/g' ./results/DESeq2/gene_logFC_pval_padj_gene_names_DOWN.txt
sed -i 's/ /\t/g' ./results/DESeq2/gene_logFC_pval_padj_gene_names_ALL.txt
rm ./results/DESeq2/gene_logFC_pval_padj_ENSGonly_UP.txt ./results/DESeq2/gene_logFC_pval_padj_ENSGonly_DOWN.txt ./results/DESeq2/gene_logFC_pval_padj_ENSGonly_ALL.txt

sed 's/\./,/g' ./results/DESeq2/gene_logFC_pval_padj_gene_names_UP.txt | sort -gk5 > ./results/DESeq2/gene_logFC_pval_padj_gene_names_padj-sorted_UP.txt
sed 's/\./,/g' ./results/DESeq2/gene_logFC_pval_padj_gene_names_DOWN.txt | sort -gk5 > ./results/DESeq2/gene_logFC_pval_padj_gene_names_padj-sorted_DOWN.txt
sed 's/\./,/g' ./results/DESeq2/gene_logFC_pval_padj_gene_names_ALL.txt | sort -gk5 > ./results/DESeq2/gene_logFC_pval_padj_gene_names_padj-sorted_ALL.txt
rm ./results/DESeq2/gene_logFC_pval_padj_gene_names_UP.txt ./results/DESeq2/gene_logFC_pval_padj_gene_names_DOWN.txt ./results/DESeq2/gene_logFC_pval_padj_gene_names_ALL.txt

tail -n +2 ./results/DESeq2/gene_logFC_pval_padj_gene_names_padj-sorted_UP.txt | cut -f1 > ./results/DESeq2/gene_names_UP.txt
tail -n +2 ./results/DESeq2/gene_logFC_pval_padj_gene_names_padj-sorted_DOWN.txt | cut -f1 > ./results/DESeq2/gene_names_DOWN.txt
