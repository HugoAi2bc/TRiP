#! /bin/bash

#####
# Save gene names, padj and logFC after DESeq2 analysis
#####


# awk '{print $1"\t"$13"\t"$14"\t"$16"\t"$17;}' ./RESULTS/DESeq2/up.txt > ./RESULTS/DESeq2/tmp.txt
# tail -n +2 ./RESULTS/DESeq2/tmp.txt > ./RESULTS/DESeq2/test_gene_logFC_pval_padj_UP.txt
#
# sed -r -e 's/:.+:CDS\t/\t/g' ./RESULTS/DESeq2/test_gene_logFC_pval_padj_UP.txt > ./RESULTS/DESeq2/test_gene_logFC_pval_padj_ENSGonly_UP.txt
#
# head -1 ./RESULTS/DESeq2/test_gene_logFC_pval_padj_UP.txt > ./RESULTS/DESeq2/test_gene_logFC_pval_padj_gene_names_UP.txt
# # FAUT TRIER POUR LE JOIN !!!
# join -2 2 -o 2.3,1.2,1.3,1.4,1.5 --nocheck-order ./RESULTS/DESeq2/test_gene_logFC_pval_padj_ENSGonly_UP.txt ./correspondancetranscript_GeneID_GeneName.sort.txt >>\
# ./RESULTS/DESeq2/test_gene_logFC_pval_padj_gene_names_UP.txt
#
# sed -i 's/ /\t/g' ./RESULTS/DESeq2/test_gene_logFC_pval_padj_gene_names_UP.txt
#
# tail -n +2 ./RESULTS/DESeq2/test_gene_logFC_pval_padj_gene_names_UP.txt | cut -f1 > ./RESULTS/DESeq2/test_gene_names_UP.txt



# Take infos of : Gene_name logFC lfcSE pval padj
awk '{print $1"\t"$13"\t"$14"\t"$16"\t"$17;}' ./RESULTS/DESeq2/up.txt > ./RESULTS/DESeq2/tmp.txt
# Sort on padj
(head -n 1 ./RESULTS/DESeq2/tmp.txt && tail -n +2 ./RESULTS/DESeq2/tmp.txt | sort -k5) > ./RESULTS/DESeq2/gene_logFC_pval_padj_UP.txt
rm ./RESULTS/DESeq2/tmp.txt
awk '{print $1"\t"$13"\t"$14"\t"$16"\t"$17;}' ./RESULTS/DESeq2/down.txt > ./RESULTS/DESeq2/tmp.txt
(head -n 1 ./RESULTS/DESeq2/tmp.txt && tail -n +2 ./RESULTS/DESeq2/tmp.txt | sort -k5) > ./RESULTS/DESeq2/gene_logFC_pval_padj_DOWN.txt
rm ./RESULTS/DESeq2/tmp.txt
awk '{print $1"\t"$13"\t"$14"\t"$16"\t"$17;}' ./RESULTS/DESeq2/complete.txt > ./RESULTS/DESeq2/tmp.txt
(head -n 1 ./RESULTS/DESeq2/tmp.txt && tail -n +2 ./RESULTS/DESeq2/tmp.txt | sort -k5) > ./RESULTS/DESeq2/gene_logFC_pval_padj_ALL.txt
rm ./RESULTS/DESeq2/tmp.txt

tail -n +2 ./RESULTS/DESeq2/gene_logFC_pval_padj_UP.txt | sed -r -e 's/:.+:CDS\t/\t/g' | sort > ./RESULTS/DESeq2/gene_logFC_pval_padj_ENSGonly_UP.txt
tail -n +2 ./RESULTS/DESeq2/gene_logFC_pval_padj_DOWN.txt | sed -r -e 's/:.+:CDS\t/\t/g' | sort > ./RESULTS/DESeq2/gene_logFC_pval_padj_ENSGonly_DOWN.txt
tail -n +2 ./RESULTS/DESeq2/gene_logFC_pval_padj_ALL.txt | sed -r -e 's/:.+:CDS\t/\t/g' | sort > ./RESULTS/DESeq2/gene_logFC_pval_padj_ENSGonly_ALL.txt

head -1 ./RESULTS/DESeq2/gene_logFC_pval_padj_UP.txt > ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_UP.txt
join -2 2 -o 2.3,1.2,1.3,1.4,1.5 ./RESULTS/DESeq2/gene_logFC_pval_padj_ENSGonly_UP.txt ./correspondancetranscript_GeneID_GeneName.sort.txt >>\
./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_UP.txt
head -1 ./RESULTS/DESeq2/gene_logFC_pval_padj_DOWN.txt > ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_DOWN.txt
join -2 2 -o 2.3,1.2,1.3,1.4,1.5 ./RESULTS/DESeq2/gene_logFC_pval_padj_ENSGonly_DOWN.txt ./correspondancetranscript_GeneID_GeneName.sort.txt >>\
./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_DOWN.txt
head -1 ./RESULTS/DESeq2/gene_logFC_pval_padj_ALL.txt > ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_ALL.txt
join -2 2 -o 2.3,1.2,1.3,1.4,1.5 ./RESULTS/DESeq2/gene_logFC_pval_padj_ENSGonly_ALL.txt ./correspondancetranscript_GeneID_GeneName.sort.txt >>\
./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_ALL.txt
rm ./RESULTS/DESeq2/gene_logFC_pval_padj_UP.txt ./RESULTS/DESeq2/gene_logFC_pval_padj_DOWN.txt ./RESULTS/DESeq2/gene_logFC_pval_padj_ALL.txt

sed -i 's/ /\t/g' ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_UP.txt
sed -i 's/ /\t/g' ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_DOWN.txt
sed -i 's/ /\t/g' ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_ALL.txt
rm ./RESULTS/DESeq2/gene_logFC_pval_padj_ENSGonly_UP.txt ./RESULTS/DESeq2/gene_logFC_pval_padj_ENSGonly_DOWN.txt ./RESULTS/DESeq2/gene_logFC_pval_padj_ENSGonly_ALL.txt

sed 's/\./,/g' ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_UP.txt | sort -gk5 > ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_padj-sorted_UP.txt
sed 's/\./,/g' ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_DOWN.txt | sort -gk5 > ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_padj-sorted_DOWN.txt
sed 's/\./,/g' ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_ALL.txt | sort -gk5 > ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_padj-sorted_ALL.txt
rm ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_UP.txt ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_DOWN.txt ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_ALL.txt

tail -n +2 ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_padj-sorted_UP.txt | cut -f1 > ./RESULTS/DESeq2/gene_names_UP.txt
tail -n +2 ./RESULTS/DESeq2/gene_logFC_pval_padj_gene_names_padj-sorted_DOWN.txt | cut -f1 > ./RESULTS/DESeq2/gene_names_DOWN.txt
