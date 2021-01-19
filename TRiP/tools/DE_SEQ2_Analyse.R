rm(list=ls())
library(DESeq2)
library(stringr)


args = commandArgs(trailingOnly=TRUE)

if (is.na(args[1])) {args[1] <- "CT"}
if (is.na(args[2])) {args[2] <- 0}
if (is.na(args[3])) {args[3] <- 0.01}

Var_log2FC <- as.numeric(args[2])
Var_padj <- as.numeric(args[3])

# Reading the data file
# expData = read.table("TEST1_RPL10.txt", header = T, row.names = 1)
#expData = read.table(snakemake@input[["counts_matrix"]], header = T, row.names = 1)
expData = read.table("./results/DESeq2/count_matrix.txt", header = T, row.names = 1)


# determination du nombre des échantillons
expData_Sorted <- expData[,sort(colnames(expData))]
Names_Col <- colnames(expData_Sorted)
Cond <- str_replace(Names_Col,regex("[[:digit:]]{1,}$"),"")



nb_col_condA <- table(factor(Cond))[[1]]
nb_col_condB <- table(factor(Cond))[[2]]



if (Cond[1] != args[1]) {
        Ord <- c( (1:nb_col_condB)+nb_col_condA , (1:nb_col_condA) )
        expData_Sorted <- expData_Sorted[, Ord]

# inversion des 2 valeurs nb_col_condA et B sans passer par une 3emme variable (lol)
        nb_col_condA <- nb_col_condA + nb_col_condB
        nb_col_condB <- nb_col_condA - nb_col_condB
        nb_col_condA <- nb_col_condA - nb_col_condB

}


# Comparaison des tailles des librairies
barplot(colSums(expData_Sorted), ylab = "Total read number",
        main = "Library sizes", col = c(rep("yellow", nb_col_condA), rep("blue", nb_col_condB)),
        names = colnames(expData_Sorted),
        cex.names = 0.6
        )

# Normalisation avec DEseq2.

conds    = factor(c(rep("CondA", nb_col_condA), rep("CondB", nb_col_condB)))

colData  = data.frame(condition = conds)

ddsObjet = DESeqDataSetFromMatrix(countData = expData_Sorted,
                                  colData   = colData, formula(~ condition))

ddsObjet = estimateSizeFactors(ddsObjet)

Size_Factors <- sizeFactors(ddsObjet)


# Données de comptage normalisées :
normCountData = counts(ddsObjet, normalized = TRUE)


# Comparaison des tailles des librairies (après normalisation)
barplot(colSums(normCountData), ylab = "Total read number",
        main = "Library sizes (after normalization)",
        col = c(rep("yellow", nb_col_condA), rep("blue", nb_col_condB)),
        names = colnames(expData_Sorted),
        cex.names = 0.6)



# La valeur +1 est ajoutée à tous les gènes, pour toutes les conditions,
# afin d'éviter les éventuelles divisions par 0.
normCountData = normCountData + 1

# Mean
CT_Mean<- data.frame(apply(normCountData[,1:(nb_col_condA)],1,mean))
Mut_Mean<- data.frame(apply(normCountData[,(nb_col_condA+1):(nb_col_condB)],1,mean))



logFC = log2(Mut_Mean / CT_Mean)
colnames(logFC) = paste("logFC", 1, sep = "_")

#---------------------------------------------
# Calcul pour chaque gène, la moyenne des logFC
# (compilation des réplicats).
# l'histogramme des valeurs obtenues.
#---------------------------------------------

meanLogFC = apply(logFC, 1, mean)

hist(meanLogFC, nclass = 100, xlab = "logFC values")

abline(v = 2, col = "red")
abline(v = -2, col = "red")


#---------------------------------------------
# L'écart type des logFC pour chaque gène
# Tracer le graphique présenté en cours (logFC en fonction de SD)
#---------------------------------------------

# sdLogFC = apply(logFC, 1, sd)
# plot(sdLogFC, meanLogFC, pch = 20)
# abline(h = 2, col = "red")
# abline(h = -2, col = "red")

#---------------------------------------------
# Calcul pour chaque gène, la valeur T, l'histogramme
# des valeurs obtenues
#---------------------------------------------

# tVal = meanLogFC / (sdLogFC/sqrt(3))
# hist(tVal, nclass = 100, xlab = "T values")


#---------------------------------------------
# Calcul pour chaque gène la moyenne des mesures
# de comptage
#---------------------------------------------

meanExp = apply(normCountData, 1, mean)

#---------------------------------------------
#  Graph des meanLogFC en fonction
# de meanExp
#---------------------------------------------


plot(meanExp[names(meanLogFC)], meanLogFC,
     pch = 20,
     log = "x",
     )
abline(h = 2, col = "red")
abline(h = -2, col = "red")


# Application du package DESeq2, pour le calcul des statistiques
ddsEstim = DESeq(ddsObjet)
resDESeq = results(ddsEstim, contrast = c("condition", "CondB", "CondA"))


# Pour obtenir une description des colonnes
mcols(resDESeq, use.names = TRUE)

# Représentation MA plot
plotMA(resDESeq)


# Pour regarder le paramètre de dispersion
dds = estimateDispersions(ddsObjet)
plotDispEsts(dds)

#---------------------------------------------
# l'histogramme des pvalues
#---------------------------------------------


hist(resDESeq[,"pvalue"], nclass = 100, xlab = "pvalue",
     main = "Adjusted p-values (DESeq2)")

#---------------------------------------------
# l'histogramme des pvalues ajustées
#---------------------------------------------


hist(resDESeq[,"padj"], nclass = 100, xlab = "padj",
     main = "Adjusted p-values (DESeq2)")
#---------------------------------------------
#  le volcano plot
#---------------------------------------------
myColors <- ifelse((resDESeq$padj < 0.01 & resDESeq$log2FoldChange >0), "green" ,
                   ifelse((resDESeq$padj < 0.01 & resDESeq$log2FoldChange <0), "red" ,
                          "black" ) )

plot(resDESeq[, "log2FoldChange"], -log(resDESeq[, "padj"]),
     pch = 20,
     col=myColors,
     xlab = "logFC", ylab = "-log(padj)")



Bruts_Norm <- cbind(expData_Sorted, normCountData )

#---------------------------------------------
# Création Tables finales , complete, up and down
#---------------------------------------------




Table_Complete <- data.frame(cbind(row.names(Bruts_Norm) , Bruts_Norm, CT_Mean,Mut_Mean , resDESeq) )
colnames(Table_Complete) = col.names = c("genes" ,
                               paste("CT", 1:nb_col_condA, sep = "_"),
                               paste("Mut", 1:nb_col_condB, sep = "_"),
                               paste("norm.CT", 1:nb_col_condA, sep = "_"),
                               paste("norm.Mut", 1:nb_col_condB, sep = "_"),
                               "CT_Mean", "Mut_Mean",
                               colnames(resDESeq))


inducedGenes = Table_Complete[which((Table_Complete[, "log2FoldChange"] > Var_log2FC) & (Table_Complete[, "padj"] < Var_padj)),]
dim(inducedGenes)

repressedGenes = Table_Complete[which((Table_Complete[, "log2FoldChange"] < -Var_log2FC) & (Table_Complete[, "padj"] < Var_padj)),]
dim(repressedGenes)

# Tri des tables par padj
inducedGenes <- inducedGenes[(order(inducedGenes[,"padj"])),]
repressedGenes <- repressedGenes[(order(repressedGenes[,"padj"])),]


# Ecriture des résultats
write.table(Table_Complete,
            file = paste0("./results/DESeq2/complete.txt"),
#            file = paste0(Cond[[nb_col_condB]],"_",Cond[[nb_col_condA+nb_col_condB]],".complete.txt"),
            quote = F, sep = "\t", row.names = F)


write.table(inducedGenes,
            file = paste0("./results/DESeq2/up.txt"),
#            file = paste0(Cond[[nb_col_condB]],"_",Cond[[nb_col_condA+nb_col_condB]],".up.txt"),
            quote = F, sep = "\t", row.names = F)

write.table(repressedGenes,
            file = paste0("./results/DESeq2/down.txt"),
#            file = paste0(Cond[[nb_col_condB]],"_",Cond[[nb_col_condA+nb_col_condB]],".down.txt"),
            quote = F, sep = "\t", row.names = F)
