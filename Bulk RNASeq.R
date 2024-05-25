### Data intake
rm(list = ls())

library(edgeR)
library(ggplot2)
library(biomaRt)
library(stringr)
library(dplyr)
library(forcats)
library (data.table)
library(tidyverse)


setwd("~/Documents/2023/R21_Malaria_study/RNAseq/Analysis/")

# Read count matrix
counts <- read.table("Export/Dmitri/Count.table.annotated.txt", sep = "\t", header = TRUE) #What does the last two things do?
ann <- read.table("Export/Dmitri/Metadata.txt", sep = "\t", header = TRUE) #What does the last two things do?

# separate count matrix into groups, eval = FALSE
g3.counts <- counts[, colnames(counts) %in% ann$Sample_ID[ann$Group == 3]]
g4.counts <- counts[, colnames(counts) %in% ann$Sample_ID[ann$Group == 4]]
g5.counts <- counts[, colnames(counts) %in% ann$Sample_ID[ann$Group == 5]]
g6.counts <- counts[, colnames(counts) %in% ann$Sample_ID[ann$Group == 6]]

ann.3 <- subset(ann, Group == 3)
ann.4 <- subset(ann, Group == 4)
ann.5 <- subset(ann, Group == 5)
ann.6 <- subset(ann, Group == 6)


#Create DGE objects and filter, eval = FALSE ######## Can you explain this?
g3 <- DGEList(counts = g3.counts, genes = counts$Gene.symbol)
g4 <- DGEList(counts = g4.counts, genes = counts$Gene.symbol)
g5 <- DGEList(counts = g5.counts, genes = counts$Gene.symbol)
g6 <- DGEList(counts = g6.counts, genes = counts$Gene.symbol)

keep.g3 <- rowSums(cpm(g3) > 1) >= floor(ncol(g3)/5) ##ceil makes it round up
keep.g4 <- rowSums(cpm(g4) > 1) >= floor(ncol(g4)/5)
keep.g5 <- rowSums(cpm(g5) > 1) >= floor(ncol(g5)/5)
keep.g6 <- rowSums(cpm(g6) > 1) >= floor(ncol(g6)/5)

g3 <- g3[keep.g3, ] # This leaves us with 11,449 genes
g4 <- g4[keep.g4, ] # This leaves us with 11,451genes
g5 <- g5[keep.g5, ] # This leaves us with 11,451 genes
g6 <- g6[keep.g6, ] # This leaves us with 11,456 genes


## DGE significance testing, eval = FALSE}

# Arrange the samples in the annotation tables in the same order as the columns in count matrices 
ann.3 <- ann.3[match(row.names(g3$samples), ann.3$Sample_ID), ]
ann.4 <- ann.4[match(row.names(g4$samples), ann.4$Sample_ID), ]
ann.5 <- ann.5[match(row.names(g5$samples), ann.5$Sample_ID), ]
ann.6 <- ann.6[match(row.names(g6$samples), ann.6$Sample_ID), ]

# Further separate counts and annotations into post-prime and post-tertiary samples
ann.3.prime <- subset(ann.3, (Time == "Baseline.Day.0" | Time == "W0.Day.1" | Time == "W1.Day.0" | Time == "W4.Day.0"))
ann.3.tert <- subset(ann.3, (Time == "W8.Day.0" | Time == "W8.Day.1" | Time == "W9.Day.0" | Time == "W12.Day.0"))

ann.4.prime <- subset(ann.4, (Time == "Baseline.Day.0" | Time == "W0.Day.1" | Time == "W1.Day.0" | Time == "W4.Day.0"))
ann.4.tert <- subset(ann.4, (Time == "W8.Day.0" | Time == "W8.Day.1" | Time == "W9.Day.0" | Time == "W12.Day.0"))

ann.5.prime <- subset(ann.5, (Time == "Baseline.Day.0" | Time == "W0.Day.1" | Time == "W1.Day.0" | Time == "W4.Day.0"))
ann.5.tert <- subset(ann.5, (Time == "W8.Day.0" | Time == "W8.Day.1" | Time == "W9.Day.0" | Time == "W12.Day.0"))

ann.6.prime <- subset(ann.6, (Time == "Baseline.Day.0" | Time == "W0.Day.1" | Time == "W1.Day.0" | Time == "W4.Day.0"))
ann.6.tert <- subset(ann.6, (Time == "W23.Day.0" | Time == "W23.Day.1" | Time == "W24.Day.0" | Time == "W27.Day.0"))

g3.prime <- g3[, row.names(g3$samples) %in% ann.3.prime$Sample_ID]
g3.tert <- g3[, row.names(g3$samples) %in% ann.3.tert$Sample_ID]

g4.prime <- g4[, row.names(g4$samples) %in% ann.4.prime$Sample_ID]
g4.tert <- g4[, row.names(g4$samples) %in% ann.4.tert$Sample_ID]

g5.prime <- g5[, row.names(g5$samples) %in% ann.5.prime$Sample_ID]
g5.tert <- g5[, row.names(g5$samples) %in% ann.5.tert$Sample_ID]

g6.prime <- g6[, row.names(g6$samples) %in% ann.6.prime$Sample_ID]
g6.tert <- g6[, row.names(g6$samples) %in% ann.6.tert$Sample_ID]

ann.3.prime$Sample_ID == row.names(g3.prime$samples) # ALL TRUE
ann.4.prime$Sample_ID == row.names(g4.prime$samples) # ALL TRUE
ann.5.prime$Sample_ID == row.names(g5.prime$samples) # ALL TRUE
ann.6.prime$Sample_ID == row.names(g6.prime$samples) # ALL TRUE

ann.3.tert$Sample_ID == row.names(g3.tert$samples) # ALL TRUE
ann.4.tert$Sample_ID == row.names(g4.tert$samples) # ALL TRUE
ann.5.tert$Sample_ID == row.names(g5.tert$samples) # ALL TRUE
ann.6.tert$Sample_ID == row.names(g6.tert$samples) # ALL TRUE


# Prepare data for making design matrices 
animal.3.prime <- ann.3.prime$Animal
animal.4.prime <- ann.4.prime$Animal
animal.5.prime <- ann.5.prime$Animal
animal.6.prime <- ann.6.prime$Animal

animal.3.tert <- ann.3.tert$Animal
animal.4.tert <- ann.4.tert$Animal
animal.5.tert <- ann.5.tert$Animal
animal.6.tert <- ann.6.tert$Animal

time.3.prime <- factor(ann.3.prime$Time, levels = c("Baseline.Day.0", "W0.Day.1", "W1.Day.0", "W4.Day.0"))
time.4.prime <- factor(ann.4.prime$Time, levels = c("Baseline.Day.0", "W0.Day.1", "W1.Day.0", "W4.Day.0"))
time.5.prime <- factor(ann.5.prime$Time, levels = c("Baseline.Day.0", "W0.Day.1", "W1.Day.0", "W4.Day.0"))
time.6.prime <- factor(ann.6.prime$Time, levels = c("Baseline.Day.0", "W0.Day.1", "W1.Day.0", "W4.Day.0"))

time.3.tert <- factor(ann.3.tert$Time, levels = c("W8.Day.0", "W8.Day.1", "W9.Day.0", "W12.Day.0"))
time.4.tert <- factor(ann.4.tert$Time, levels = c("W8.Day.0", "W8.Day.1", "W9.Day.0", "W12.Day.0"))
time.5.tert <- factor(ann.5.tert$Time, levels = c("W8.Day.0", "W8.Day.1", "W9.Day.0", "W12.Day.0"))
time.6.tert <- factor(ann.6.tert$Time, levels = c("W23.Day.0", "W23.Day.1", "W24.Day.0", "W27.Day.0"))

# Define the design matrices for comparing everything to the respective baselines
design.g3.prime <- model.matrix(~animal.3.prime + time.3.prime, data=g3.prime$samples)
design.g4.prime <- model.matrix(~animal.4.prime + time.4.prime, data=g4.prime$samples)
design.g5.prime <- model.matrix(~animal.5.prime + time.5.prime, data=g5.prime$samples)
design.g6.prime <- model.matrix(~animal.6.prime + time.6.prime, data=g6.prime$samples)

design.g3.tert <- model.matrix(~animal.3.tert + time.3.tert, data=g3.tert$samples)
design.g4.tert <- model.matrix(~animal.4.tert + time.4.tert, data=g4.tert$samples)
design.g5.tert <- model.matrix(~animal.5.tert + time.5.tert, data=g5.tert$samples)
design.g6.tert <- model.matrix(~animal.6.tert + time.6.tert, data=g6.tert$samples)

# Estimate dispersions
g3.prime <- estimateDisp(g3.prime, design.g3.prime)
g4.prime <- estimateDisp(g4.prime, design.g4.prime)
g5.prime <- estimateDisp(g5.prime, design.g5.prime)
g6.prime <- estimateDisp(g6.prime, design.g6.prime)

g3.tert <- estimateDisp(g3.tert, design.g3.tert)
g4.tert <- estimateDisp(g4.tert, design.g4.tert)
g5.tert <- estimateDisp(g5.tert, design.g5.tert)
g6.tert <- estimateDisp(g6.tert, design.g6.tert)

# Caclulate normalization factors
g3.prime <- calcNormFactors(g3.prime) #What is this?
g4.prime <- calcNormFactors(g4.prime)
g5.prime <- calcNormFactors(g5.prime)
g6.prime <- calcNormFactors(g6.prime)

g3.tert <- calcNormFactors(g3.tert)
g4.tert <- calcNormFactors(g4.tert)
g5.tert <- calcNormFactors(g5.tert)
g6.tert <- calcNormFactors(g6.tert)


# Fit models
fit.g3.prime <- glmQLFit(g3.prime, design.g3.prime)
fit.g4.prime <- glmQLFit(g4.prime, design.g4.prime)
fit.g5.prime <- glmQLFit(g5.prime, design.g5.prime)
fit.g6.prime <- glmQLFit(g6.prime, design.g6.prime)

fit.g3.tert <- glmQLFit(g3.tert, design.g3.tert)
fit.g4.tert <- glmQLFit(g4.tert, design.g4.tert)
fit.g5.tert <- glmQLFit(g5.tert, design.g5.tert)
fit.g6.tert <- glmQLFit(g6.tert, design.g6.tert)

# Calculate significance
lrt_g3.prime.d1 <- glmQLFTest(fit.g3.prime, coef = 7) #What are these numbers? This comes from the deisgn matrix ##colnames(design.g3.prime)
lrt_g3.prime.d7 <- glmQLFTest(fit.g3.prime, coef = 8)
lrt_g3.prime.d28 <- glmQLFTest(fit.g3.prime, coef = 9)

lrt_g4.prime.d1 <- glmQLFTest(fit.g4.prime, coef = 7)
lrt_g4.prime.d7 <- glmQLFTest(fit.g4.prime, coef = 8)
lrt_g4.prime.d28 <- glmQLFTest(fit.g4.prime, coef = 9)

lrt_g5.prime.d1 <- glmQLFTest(fit.g5.prime, coef = 11)
lrt_g5.prime.d7 <- glmQLFTest(fit.g5.prime, coef = 12)
lrt_g5.prime.d28 <- glmQLFTest(fit.g5.prime, coef = 13)

lrt_g6.prime.d1 <- glmQLFTest(fit.g6.prime, coef = 11)
lrt_g6.prime.d7 <- glmQLFTest(fit.g6.prime, coef = 12)
lrt_g6.prime.d28 <- glmQLFTest(fit.g6.prime, coef = 13)

lrt_g3.tert.d1 <- glmQLFTest(fit.g3.tert, coef = 7)
lrt_g3.tert.d7 <- glmQLFTest(fit.g3.tert, coef = 8)
lrt_g3.tert.d28 <- glmQLFTest(fit.g3.tert, coef = 9)

lrt_g4.tert.d1 <- glmQLFTest(fit.g4.tert, coef = 7)
lrt_g4.tert.d7 <- glmQLFTest(fit.g4.tert, coef = 8)
lrt_g4.tert.d28 <- glmQLFTest(fit.g4.tert, coef = 9)

lrt_g5.tert.d1 <- glmQLFTest(fit.g5.tert, coef = 11)
lrt_g5.tert.d7 <- glmQLFTest(fit.g5.tert, coef = 12)
lrt_g5.tert.d28 <- glmQLFTest(fit.g5.tert, coef = 13)

lrt_g6.tert.d1 <- glmQLFTest(fit.g6.tert, coef = 11)
lrt_g6.tert.d7 <- glmQLFTest(fit.g6.tert, coef = 12)
lrt_g6.tert.d28 <- glmQLFTest(fit.g6.tert, coef = 13)


g3.prime.d1.degs <- as.data.frame(topTags(lrt_g3.prime.d1, n=Inf))
g3.prime.d7.degs <- as.data.frame(topTags(lrt_g3.prime.d7, n=Inf))
g3.prime.d28.degs <- as.data.frame(topTags(lrt_g3.prime.d28, n=Inf))

g4.prime.d1.degs <- as.data.frame(topTags(lrt_g4.prime.d1, n=Inf))
g4.prime.d7.degs <- as.data.frame(topTags(lrt_g4.prime.d7, n=Inf))
g4.prime.d28.degs <- as.data.frame(topTags(lrt_g4.prime.d28, n=Inf))

g5.prime.d1.degs <- as.data.frame(topTags(lrt_g5.prime.d1, n=Inf))
g5.prime.d7.degs <- as.data.frame(topTags(lrt_g5.prime.d7, n=Inf))
g5.prime.d28.degs <- as.data.frame(topTags(lrt_g5.prime.d28, n=Inf))

g6.prime.d1.degs <- as.data.frame(topTags(lrt_g6.prime.d1, n=Inf))
g6.prime.d7.degs <- as.data.frame(topTags(lrt_g6.prime.d7, n=Inf))
g6.prime.d28.degs <- as.data.frame(topTags(lrt_g6.prime.d28, n=Inf))

g3.tert.d1.degs <- as.data.frame(topTags(lrt_g3.tert.d1, n=Inf))
g3.tert.d7.degs <- as.data.frame(topTags(lrt_g3.tert.d7, n=Inf))
g3.tert.d28.degs <- as.data.frame(topTags(lrt_g3.tert.d28, n=Inf))

g4.tert.d1.degs <- as.data.frame(topTags(lrt_g4.tert.d1, n=Inf))
g4.tert.d7.degs <- as.data.frame(topTags(lrt_g4.tert.d7, n=Inf))
g4.tert.d28.degs <- as.data.frame(topTags(lrt_g4.tert.d28, n=Inf))

g5.tert.d1.degs <- as.data.frame(topTags(lrt_g5.tert.d1, n=Inf))
g5.tert.d7.degs <- as.data.frame(topTags(lrt_g5.tert.d7, n=Inf))
g5.tert.d28.degs <- as.data.frame(topTags(lrt_g5.tert.d28, n=Inf))

g6.tert.d1.degs <- as.data.frame(topTags(lrt_g6.tert.d1, n=Inf))
g6.tert.d7.degs <- as.data.frame(topTags(lrt_g6.tert.d7, n=Inf))
g6.tert.d28.degs <- as.data.frame(topTags(lrt_g6.tert.d28, n=Inf))

###Significant DEGs
########Group 3
g3.prime.d1.sigpos <- g3.prime.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g3.prime.d1.signeg <- g3.prime.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g3.prime.d1.allsig <- rbind (g3.prime.d1.sigpos, g3.prime.d1.signeg)

g3.tert.d1.sigpos <- g3.tert.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g3.tert.d1.signeg <- g3.tert.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g3.tert.d1.allsig <- rbind (g3.tert.d1.sigpos, g3.tert.d1.signeg)

g3.prime.d7.sigpos <- g3.prime.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g3.prime.d7.signeg <- g3.prime.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g3.prime.d7.allsig <- rbind (g3.prime.d7.sigpos, g3.prime.d7.signeg)

g3.tert.d7.sigpos <- g3.tert.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g3.tert.d7.signeg <- g3.tert.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g3.tert.d7.allsig <- rbind (g3.tert.d7.sigpos, g3.tert.d7.signeg)

g3.prime.d28.sigpos <- g3.prime.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g3.prime.d28.signeg <- g3.prime.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g3.prime.d28.allsig <- rbind (g3.prime.d28.sigpos, g3.prime.d28.signeg)

g3.tert.d28.sigpos <- g3.tert.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g3.tert.d28.signeg <- g3.tert.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g3.tert.d28.allsig <- rbind (g3.tert.d28.sigpos, g3.tert.d28.signeg)


########Group 4
g4.prime.d1.sigpos <- g4.prime.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g4.prime.d1.signeg <- g4.prime.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g4.prime.d1.allsig <- rbind (g4.prime.d1.sigpos, g4.prime.d1.signeg)

g4.tert.d1.sigpos <- g4.tert.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g4.tert.d1.signeg <- g4.tert.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g4.tert.d1.allsig <- rbind (g4.tert.d1.sigpos, g4.tert.d1.signeg)

g4.prime.d7.sigpos <- g4.prime.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g4.prime.d7.signeg <- g4.prime.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g4.prime.d7.allsig <- rbind (g4.prime.d7.sigpos, g4.prime.d7.signeg)

g4.tert.d7.sigpos <- g4.tert.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g4.tert.d7.signeg <- g4.tert.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g4.tert.d7.allsig <- rbind (g4.tert.d7.sigpos, g4.tert.d7.signeg)

g4.prime.d28.sigpos <- g4.prime.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g4.prime.d28.signeg <- g4.prime.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g4.prime.d28.allsig <- rbind (g4.prime.d28.sigpos, g4.prime.d28.signeg)

g4.tert.d28.sigpos <- g4.tert.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g4.tert.d28.signeg <- g4.tert.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g4.tert.d28.allsig <- rbind (g4.tert.d28.sigpos, g4.tert.d28.signeg)


########Group 5
g5.prime.d1.sigpos <- g5.prime.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g5.prime.d1.signeg <- g5.prime.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g5.prime.d1.allsig <- rbind (g5.prime.d1.sigpos, g5.prime.d1.signeg)

g5.tert.d1.sigpos <- g5.tert.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g5.tert.d1.signeg <- g5.tert.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g5.tert.d1.allsig <- rbind (g5.tert.d1.sigpos, g5.tert.d1.signeg)

g5.prime.d7.sigpos <- g5.prime.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g5.prime.d7.signeg <- g5.prime.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g5.prime.d7.allsig <- rbind (g5.prime.d7.sigpos, g5.prime.d7.signeg)

g5.tert.d7.sigpos <- g5.tert.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g5.tert.d7.signeg <- g5.tert.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g5.tert.d7.allsig <- rbind (g5.tert.d7.sigpos, g5.tert.d7.signeg)

g5.prime.d28.sigpos <- g5.prime.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g5.prime.d28.signeg <- g5.prime.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g5.prime.d28.allsig <- rbind (g5.prime.d28.sigpos, g5.prime.d28.signeg)

g5.tert.d28.sigpos <- g5.tert.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g5.tert.d28.signeg <- g5.tert.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g5.tert.d28.allsig <- rbind (g5.tert.d28.sigpos, g5.tert.d28.signeg)


########Group 6
g6.prime.d1.sigpos <- g6.prime.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g6.prime.d1.signeg <- g6.prime.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g6.prime.d1.allsig <- rbind (g6.prime.d1.sigpos, g6.prime.d1.signeg)

g6.tert.d1.sigpos <- g6.tert.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g6.tert.d1.signeg <- g6.tert.d1.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g6.tert.d1.allsig <- rbind (g6.tert.d1.sigpos, g6.tert.d1.signeg)

g6.prime.d7.sigpos <- g6.prime.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g6.prime.d7.signeg <- g6.prime.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g6.prime.d7.allsig <- rbind (g6.prime.d7.sigpos, g6.prime.d7.signeg)

g6.tert.d7.sigpos <- g6.tert.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g6.tert.d7.signeg <- g6.tert.d7.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g6.tert.d7.allsig <- rbind (g6.tert.d7.sigpos, g6.tert.d7.signeg)

g6.prime.d28.sigpos <- g6.prime.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g6.prime.d28.signeg <- g6.prime.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g6.prime.d28.allsig <- rbind (g6.prime.d28.sigpos, g6.prime.d28.signeg)

g6.tert.d28.sigpos <- g6.tert.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC > 1)
g6.tert.d28.signeg <- g6.tert.d28.degs %>% 
  subset (!FDR > 0.05) %>% 
  subset (logFC < -1)
g6.tert.d28.allsig <- rbind (g6.tert.d28.sigpos, g6.tert.d28.signeg)

####### Binding all signficant ones
A01 <- nrow (g3.prime.d1.sigpos)
A02 <- nrow (g3.prime.d1.signeg)
A03 <- nrow (g3.prime.d7.sigpos)
A04 <- nrow (g3.prime.d7.signeg)
A05 <- nrow (g3.prime.d28.sigpos)
A06 <- nrow (g3.prime.d28.signeg)
A07 <- nrow (g3.tert.d1.sigpos)
A08 <- nrow (g3.tert.d1.signeg)
A09 <- nrow (g3.tert.d7.sigpos)
A10 <- nrow (g3.tert.d7.signeg)
A11 <- nrow (g3.tert.d28.sigpos)
A12 <- nrow (g3.tert.d28.signeg)

B01 <- nrow (g4.prime.d1.sigpos)
B02 <- nrow (g4.prime.d1.signeg)
B03 <- nrow (g4.prime.d7.sigpos)
B04 <- nrow (g4.prime.d7.signeg)
B05 <- nrow (g4.prime.d28.sigpos)
B06 <- nrow (g4.prime.d28.signeg)
B07 <- nrow (g4.tert.d1.sigpos)
B08 <- nrow (g4.tert.d1.signeg)
B09 <- nrow (g4.tert.d7.sigpos)
B10 <- nrow (g4.tert.d7.signeg)
B11 <- nrow (g4.tert.d28.sigpos)
B12 <- nrow (g4.tert.d28.signeg)

C01 <- nrow (g5.prime.d1.sigpos)
C02 <- nrow (g5.prime.d1.signeg)
C03 <- nrow (g5.prime.d7.sigpos)
C04 <- nrow (g5.prime.d7.signeg)
C05 <- nrow (g5.prime.d28.sigpos)
C06 <- nrow (g5.prime.d28.signeg)
C07 <- nrow (g5.tert.d1.sigpos)
C08 <- nrow (g5.tert.d1.signeg)
C09 <- nrow (g5.tert.d7.sigpos)
C10 <- nrow (g5.tert.d7.signeg)
C11 <- nrow (g5.tert.d28.sigpos)
C12 <- nrow (g5.tert.d28.signeg)

D01 <- nrow (g6.prime.d1.sigpos)
D02 <- nrow (g6.prime.d1.signeg)
D03 <- nrow (g6.prime.d7.sigpos)
D04 <- nrow (g6.prime.d7.signeg)
D05 <- nrow (g6.prime.d28.sigpos)
D06 <- nrow (g6.prime.d28.signeg)
D07 <- nrow (g6.tert.d1.sigpos)
D08 <- nrow (g6.tert.d1.signeg)
D09 <- nrow (g6.tert.d7.sigpos)
D10 <- nrow (g6.tert.d7.signeg)
D11 <- nrow (g6.tert.d28.sigpos)
D12 <- nrow (g6.tert.d28.signeg)

all_sig_DEG_numbers <- rbind (A01, A02, A03, A04, A05, A06, A07, A08, A09, A10, A11, A12,
                              B01, B02, B03, B04, B05, B06, B07, B08, B09, B10, B11, B12,
                              C01, C02, C03, C04, C05, C06, C07, C08, C09, C10, C11, C12,
                              D01, D02, D03, D04, D05, D06, D07, D08, D09, D10, D11, D12)

view (all_sig_DEG_numbers)

# Output the results
write.table(g3.prime.d1.degs, file = "Data_output/Group.3.prime.Day1.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g3.prime.d7.degs, file = "Data_output/Group.3.prime.Day7.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g3.prime.d28.degs, file = "Data_output/Group.3.prime.Day28.DEGs.txt", sep = "\t", row.names = FALSE)

write.table(g4.prime.d1.degs, file = "Data_output/Group.4.prime.Day1.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g4.prime.d7.degs, file = "Data_output/Group.4.prime.Day7.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g4.prime.d28.degs, file = "Data_output/Group.4.prime.Day28.DEGs.txt", sep = "\t", row.names = FALSE)

write.table(g5.prime.d1.degs, file = "Data_output/Group.5.prime.Day1.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g5.prime.d7.degs, file = "Data_output/Group.5.prime.Day7.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g5.prime.d28.degs, file = "Data_output/Group.5.prime.Day28.DEGs.txt", sep = "\t", row.names = FALSE)

write.table(g6.prime.d1.degs, file = "Data_output/Group.6.prime.Day1.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g6.prime.d7.degs, file = "Data_output/Group.6.prime.Day7.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g6.prime.d28.degs, file = "Data_output/Group.6.prime.Day28.DEGs.txt", sep = "\t", row.names = FALSE)

write.table(g3.tert.d1.degs, file = "Data_output/Group.3.tert.Day1.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g3.tert.d7.degs, file = "Data_output/Group.3.tert.Day7.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g3.tert.d28.degs, file = "Data_output/Group.3.tert.Day28.DEGs.txt", sep = "\t", row.names = FALSE)

write.table(g4.tert.d1.degs, file = "Data_output/Group.4.tert.Day1.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g4.tert.d7.degs, file = "Data_output/Group.4.tert.Day7.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g4.tert.d28.degs, file = "Data_output/Group.4.tert.Day28.DEGs.txt", sep = "\t", row.names = FALSE)

write.table(g5.tert.d1.degs, file = "Data_output/Group.5.tert.Day1.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g5.tert.d7.degs, file = "Data_output/Group.5.tert.Day7.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g5.tert.d28.degs, file = "Data_output/Group.5.tert.Day28.DEGs.txt", sep = "\t", row.names = FALSE)

write.table(g6.tert.d1.degs, file = "Data_output/Group.6.tert.Day1.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g6.tert.d7.degs, file = "Data_output/Group.6.tert.Day7.DEGs.txt", sep = "\t", row.names = FALSE)
write.table(g6.tert.d28.degs, file = "Data_output/Group.6.tert.Day28.DEGs.txt", sep = "\t", row.names = FALSE)


## make volcano plots, eval = TRUE

l <- list(list(g3.prime.d1.degs, "Group 3, Day 1 after primary"), 
          list(g3.prime.d7.degs, "Group 3, Day 7 after primary"), 
          list(g3.prime.d28.degs, "Group 3, Day 28 after primary"), 
          list(g4.prime.d1.degs, "Group 4, Day 1 after primary"), 
          list(g4.prime.d7.degs, "Group 4, Day 7 after primary"), 
          list(g4.prime.d28.degs, "Group 4, Day 28 after primary"), 
          list(g5.prime.d1.degs, "Group 5, Day 1 after primary"), 
          list(g5.prime.d7.degs, "Group 5, Day 7 after primary"), 
          list(g5.prime.d28.degs, "Group 5, Day 28 after primary"), 
          list(g6.prime.d1.degs, "Group 6, Day 1 after primary"), 
          list(g6.prime.d7.degs, "Group 6, Day 7 after primary"), 
          list(g6.prime.d28.degs, "Group 6, Day 28 after primary"), 
          list(g3.tert.d1.degs, "Group 3, Day 1 after tertiary"), 
          list(g3.tert.d7.degs, "Group 3, Day 7 after tertiary"), 
          list(g3.tert.d28.degs, "Group 3, Day 28 after tertiary"), 
          list(g4.tert.d1.degs, "Group 4, Day 1 after tertiary"), 
          list(g4.tert.d7.degs, "Group 4, Day 7 after tertiary"), 
          list(g4.tert.d28.degs, "Group 4, Day 28 after tertiary"), 
          list(g5.tert.d1.degs, "Group 5, Day 1 after tertiary"), 
          list(g5.tert.d7.degs, "Group 5, Day 7 after tertiary"), 
          list(g5.tert.d28.degs, "Group 5, Day 28 after tertiary"), 
          list(g6.tert.d1.degs, "Group 6, Day 1 after tertiary"), 
          list(g6.tert.d7.degs, "Group 6, Day 7 after tertiary"), 
          list(g6.tert.d28.degs, "Group 6, Day 28 after tertiary")
)

plott <- function(data, name) {
  data$log10.FDR <- -1*log10(data$FDR)
  data$color <- "grey"
  data$color[data$logFC > 0.5 & data$FDR < 0.05] <- "red"
  data$color[data$logFC < -0.5 & data$FDR < 0.05] <- "blue"
  data$alpha <- 0.25
  data$alpha[(data$logFC > 0.5 | data$logFC < -0.5) & data$FDR < 0.05] <- 1
  
  ggplot(data, aes(x = logFC, y = log10.FDR)) + 
    geom_point(alpha = data$alpha, color = data$color) + 
    geom_hline(yintercept = -1*log10(0.05), color = "red", linetype = "dashed") + 
    geom_vline(xintercept = -0.5, linetype = "dotted") + 
    geom_vline(xintercept = 0.5, linetype = "dotted") + 
    theme_bw() +
    labs(title = name)
  ggsave(paste0("Data_output/Plots/Volcano.plot.", name, ".pdf"))
  
}

for (i in 1:length(l)) {
  data <- l[[i]][[1]]
  name <- l[[i]][[2]]
  plott(data, name)
}



####GSEA
#Run GSEA on mean log FC for all groups and time points
# GSEA, eval = FALSE}

gsea <- function(x, s) {
  
  # Run GSEA on each sample
  setwd("~/Documents/2023/R21_Malaria_study/RNAseq/Analysis/Data_output/GSEA.on.logFC/")
  system("mkdir rnk")
  system("mkdir GSEA.results")
  system("mkdir GSEA.settings")
  system("mkdir GSEA.summary")
  
  # Write out the .rnk file
  x <- subset(x, x$genes != "")
  write.table(x, file = paste0("rnk/", s, ".rnk"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  
  # Create the settings file
  
  c <- list(c("gmx", "/Users/prabhusa/Documents/2023/R21_Malaria_study/RNAseq/Analysis/Raw_data/BTM_for_GSEA_20131008.gmt"), 
            c('scoring_scheme', 'weighted'), 
            c('nperm',	'5000'), 
            c('collapse',	'No_Collapse'), 
            c('norm',	'meandiv'), 
            c('mode', 	'Max_probe'), 
            c('rnk',	paste0(getwd(), "/rnk/", s, ".rnk")), 
            c('gui',	'false'), 
            c('rpt_label',	s), 
            c('set_min',	'10'),
            c('set_max',	'500'),
            c('plot_top_x',	'50'),
            c('include_only_symbols',	'true'),
            c('make_sets',	'true'), 
            c('rnd_seed',	'timestamp'), 
            c('zip_report',	'false'), 
            c('out',	paste0(getwd(), "/GSEA.results/"))
  )
  
  system(paste0("rm ", paste0(getwd(), "/GSEA.settings/params.txt"))) # Remove the previous parameters file, if exists
  
  lapply(c, write, paste0(getwd(), "/GSEA.settings/params_", s, ".txt"), append=TRUE, ncolumns=2, sep="\t")
  
  java.comm <- paste0("java -cp /Users/prabhusa/Documents/2023/R21_Malaria_study/RNAseq/Analysis/Raw_data/gsea-3.0.jar -Xmx5000m xtools.gsea.GseaPreranked -param_file /Users/prabhusa/Documents/2023/R21_Malaria_study/RNAseq/Analysis/Data_output/GSEA.on.logFC/GSEA.settings/params_", s, ".txt")
  
  system(java.comm)
  
  setwd("../")
  
}

l <- list(list(lrt_g3.prime.d1, "G3.prime.D1"), 
          list(lrt_g3.prime.d7, "G3.prime.D7"), 
          list(lrt_g3.prime.d28, "G3.prime.D28"), 
          list(lrt_g3.tert.d1, "G3.tert.D1"), 
          list(lrt_g3.tert.d7, "G3.tert.D7"), 
          list(lrt_g3.tert.d28, "G3.tert.D28"), 
          list(lrt_g4.prime.d1, "G4.prime.D1"), 
          list(lrt_g4.prime.d7, "G4.prime.D7"), 
          list(lrt_g4.prime.d28, "G4.prime.D28"), 
          list(lrt_g4.tert.d1, "G4.tert.D1"), 
          list(lrt_g4.tert.d7, "G4.tert.D7"), 
          list(lrt_g4.tert.d28, "G4.tert.D28"), 
          list(lrt_g5.prime.d1, "G5.prime.D1"), 
          list(lrt_g5.prime.d7, "G5.prime.D7"), 
          list(lrt_g5.prime.d28, "G5.prime.D28"), 
          list(lrt_g5.tert.d1, "G5.tert.D1"), 
          list(lrt_g5.tert.d7, "G5.tert.D7"), 
          list(lrt_g5.tert.d28, "G5.tert.D28"), 
          list(lrt_g6.prime.d1, "G6.prime.D1"), 
          list(lrt_g6.prime.d7, "G6.prime.D7"), 
          list(lrt_g6.prime.d28, "G6.prime.D28"), 
          list(lrt_g6.tert.d1, "G6.tert.D1"), 
          list(lrt_g6.tert.d7, "G6.tert.D7"), 
          list(lrt_g6.tert.d28, "G6.tert.D28")
)

for (i in 1:length(l)) {
  x <- as.data.frame(topTags(l[[i]][[1]], n = Inf))[, c(1,2)]
  s <- l[[i]][[2]]
  
  gsea(x, s)
  
}

# Read in BTM annotations
btm <- read.table('~/Documents/2023/R21_Malaria_study/RNAseq/Analysis/Raw_data/BTMs.txt', header=TRUE, sep='\t')
btm <- within(btm, NAME <- toupper(NAME))

# Define function; modify the path depending on the day
all_nes <- function(dirlist){
  
  for (i in 1:length(dirlist)) {
    
    dirname <- dirlist[i]
    
    setwd(paste("~/Documents/2023/R21_Malaria_study/RNAseq/Analysis/Data_output/GSEA.on.logFC/GSEA.results/", dirname, "/", sep=""))
    
    sample_name <- unlist(strsplit(dirname, "\\.Gsea"))[[1]]
    
    excel_pos <- list.files(path = '.', pattern='^gsea_report_for_na_pos.+.xls')
    excel_neg <- list.files(path = '.', pattern='^gsea_report_for_na_neg.+.xls')
    
    pos <- read.table(excel_pos, sep='\t', header = T)
    neg <- read.table(excel_neg, sep='\t', header = T)
    pos_neg = rbind(pos, neg)
    
    # Collect all NES regardless of significance
    if (i==1) {
      all.nes <- merge(btm, pos_neg[,c(1,6)], by='NAME')
      colnames(all.nes)[i+2] <- sample_name
    } else {
      all.nes <- merge(all.nes, pos_neg[,c(1,6)], by='NAME')
      colnames(all.nes)[i+2] <- sample_name
    }
    
    # Collect only significant NES, replacing non-significant with zeros
    if (i==1) {
      sig.nes <- merge(btm, pos_neg[,c(1,6,8)], by='NAME')
      sig.nes <- within(sig.nes, NES[FDR.q.val > 0.05] <- 0)
      sig.nes$FDR.q.val <- NULL
      colnames(sig.nes)[i+2] <- sample_name
    } else {
      sig.nes <- merge(sig.nes, pos_neg[,c(1,6,8)], by='NAME')
      sig.nes <- within(sig.nes, NES[FDR.q.val > 0.05] <- 0)
      sig.nes$FDR.q.val <- NULL
      colnames(sig.nes)[i+2] <- sample_name
    }
    
    # Make a single table with a single columns of abs(NES) with annotations
    # Actually, 2 tables for post-pre-vaccination comparisons and baseline comparisons
    
    pos_neg$sample <- sample_name
    pos_neg$ABS.NES <- abs(pos_neg$NES)
    pos_neg$Significant <- ifelse((pos_neg$FDR.q.val <= 0.05), "Yes", "No")
    pos_neg$Direction <- ifelse((pos_neg$NES > 0), "Positive", "Negative")
    pos_neg <- merge(btm, pos_neg, by="NAME")
    pos_neg <- pos_neg[,c(1,2,14,15,16,17)]
    
    if (i==1) { 
      vacc <- pos_neg
    } else {
      vacc <- rbind(vacc, pos_neg)
    }
    
  }
  
  result <- list(all.nes, sig.nes, vacc)
  
  return(result)
}


e <- list.dirs(path = "~/Documents/2023/R21_Malaria_study/RNAseq/Analysis/Data_output/GSEA.on.logFC/GSEA.results/", recursive = FALSE, full.names = FALSE)

y <- all_nes(e)


varname_str <- "GSEA.all.NES"
varname <- varname_str
assign(varname, y[[1]])
p <- y[[1]]
write.table(p, paste("~/Documents/2023/R21_Malaria_study/RNAseq/Analysis/Data_output/GSEA.on.logFC/GSEA.summary/", varname_str, ".txt", sep=""), col.names=NA, sep='\t')

varname_str <- "GSEA.sig.NES"
varname <- varname_str
assign(varname, y[[2]])
p <- y[[2]]
write.table(p, paste("~/Documents/2023/R21_Malaria_study/RNAseq/Analysis/Data_output/GSEA.on.logFC/GSEA.summary/", varname_str, ".txt", sep=""), col.names=NA, sep='\t')

varname_str <- "GSEA.vacc"
varname <- varname_str
assign(varname, y[[3]])
p <- y[[3]]
write.table(p, paste("~/Documents/2023/R21_Malaria_study/RNAseq/Analysis/Data_output/GSEA.on.logFC/GSEA.summary/", varname_str, ".txt", sep=""), col.names=NA, sep='\t')

# Plotting

# Remove TBA modules
GSEA.vacc <- subset(GSEA.vacc, SUBGROUP != "TBA")

# Remove modules with insignificant enrichment across all samples
unq.modules <- unique(GSEA.vacc$NAME)
for (i in 1:length(unq.modules)) {
  my.module <- unq.modules[i]
  my.vacc <- subset(GSEA.vacc, NAME == my.module)
  if (all(my.vacc$Significant == "No")) {
    GSEA.vacc <- subset(GSEA.vacc, NAME != my.module)
  } else {
    next 
  }
}

fill <- c("Negative"= "steelblue", "Positive" = "orangered2")

GSEA.vacc$SUBGROUP <- factor(GSEA.vacc$SUBGROUP, levels=c("INFLAMMATORY/TLR/CHEMOKINES", "INTERFERON/ANTIVIRAL SENSING", 
                                                          "DC ACTIVATION", "ANTIGEN PRESENTATION", 
                                                          "MONOCYTES", "NEUTROPHILS", "NK CELLS",
                                                          "B CELLS", "PLASMA CELLS", "T CELLS", "CELL CYCLE", 
                                                          "TRANSCRIPTION AND TRANSLATION", "ECM AND MIGRATION",
                                                          "MITOCHONDRIAL/ENERGY METABOLISM", "SIGNAL TRANSDUCTION", "PLATELETS"))

GSEA.vacc$sample <- factor(GSEA.vacc$sample, levels = c("G3.prime.D1", 
                                                        "G3.prime.D7", 
                                                        "G3.prime.D28", 
                                                        "G3.tert.D1", 
                                                        "G3.tert.D7", 
                                                        "G3.tert.D28", 
                                                        "G4.prime.D1", 
                                                        "G4.prime.D7", 
                                                        "G4.prime.D28", 
                                                        "G4.tert.D1", 
                                                        "G4.tert.D7", 
                                                        "G4.tert.D28",
                                                        "G5.prime.D1", 
                                                        "G5.prime.D7", 
                                                        "G5.prime.D28", 
                                                        "G5.tert.D1", 
                                                        "G5.tert.D7", 
                                                        "G5.tert.D28",
                                                        "G6.prime.D1", 
                                                        "G6.prime.D7", 
                                                        "G6.prime.D28", 
                                                        "G6.tert.D1", 
                                                        "G6.tert.D7", 
                                                        "G6.tert.D28"
)
)

GSEA.vacc <- GSEA.vacc[order(GSEA.vacc$SUBGROUP), ]
GSEA.vacc$number <- c(1:nrow(GSEA.vacc))

GSEA.vacc %>% 
  mutate(NAME = fct_reorder(NAME, -number)) %>%
  ggplot( aes(x = sample, y = NAME, label = NAME)) +
  geom_point(aes(size = ABS.NES, colour = Direction, alpha=Significant)) + 
  scale_alpha_discrete(range=c(0.15, 0.6)) + 
  scale_radius(trans="exp", range = c(1,15)) + 
  scale_colour_manual(values=fill) + 
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x=element_blank(), axis.title.y=element_blank()) + 
  ggtitle(label="GSEA on logFC") +
  theme(panel.background = element_blank())
ggsave("Data_output/GSEA.on.logFC/GSEA.summary/Plots/All.modules.pdf", width = 15, height = 12)


# And plot by BTM functional cluster (innate, adaptive, other)

GSEA.vacc.innate <- subset(GSEA.vacc, SUBGROUP == "ANTIGEN PRESENTATION" | SUBGROUP == "DC ACTIVATION" | SUBGROUP == "INFLAMMATORY/TLR/CHEMOKINES" | SUBGROUP == "INTERFERON/ANTIVIRAL SENSING" | SUBGROUP == "MONOCYTES" | SUBGROUP == "NEUTROPHILS" | SUBGROUP == "NK CELLS")

GSEA.vacc.innate <- GSEA.vacc.innate[order(GSEA.vacc.innate$SUBGROUP), ]
GSEA.vacc.innate$number <- c(1:nrow(GSEA.vacc.innate))

GSEA.vacc.innate %>% 
  mutate(NAME = fct_reorder(NAME, -number)) %>%
  ggplot( aes(x = sample, y = NAME, label = NAME)) +
  geom_point(aes(size = ABS.NES, colour = Direction, alpha=Significant)) + 
  scale_alpha_discrete(range=c(0.15, 0.6)) + 
  scale_radius(trans="exp", range = c(1,15)) + 
  scale_colour_manual(values=fill) + 
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x=element_blank(), axis.title.y=element_blank()) + 
  ggtitle(label="GSEA on logFC, innate modules") +
  theme(panel.background = element_blank())
ggsave("Output/GSEA.on.logFC/GSEA.summary/Plots/Innate.modules.pdf", width = 15, height = 12)


GSEA.vacc.adapt <- subset(GSEA.vacc, SUBGROUP == "B CELLS" | SUBGROUP == "PLASMA CELLS" | SUBGROUP == "T CELLS" | SUBGROUP == "CELL CYCLE")

GSEA.vacc.adapt <- GSEA.vacc.adapt[order(GSEA.vacc.adapt$SUBGROUP), ]
GSEA.vacc.adapt$number <- c(1:nrow(GSEA.vacc.adapt))

GSEA.vacc.adapt %>% 
  mutate(NAME = fct_reorder(NAME, -number)) %>%
  ggplot( aes(x = sample, y = NAME, label = NAME)) +
  geom_point(aes(size = ABS.NES, colour = Direction, alpha=Significant)) + 
  scale_alpha_discrete(range=c(0.15, 0.6)) + 
  scale_radius(trans="exp", range = c(1,15)) + 
  scale_colour_manual(values=fill) + 
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x=element_blank(), axis.title.y=element_blank()) + 
  ggtitle(label="GSEA on logFC, adaptive modules") +
  theme(panel.background = element_blank())
ggsave("Output/GSEA.on.logFC/GSEA.summary/Plots/Adaptive.modules.pdf", width = 15, height = 12)


GSEA.vacc.other <- subset(GSEA.vacc, SUBGROUP == "TRANSCRIPTION AND TRANSLATION" | SUBGROUP == "SIGNAL TRANSDUCTION" | SUBGROUP == "MITOCHONDRIAL/ENERGY METABOLISM" | SUBGROUP == "ECM AND MIGRATION" | SUBGROUP == "PLATELETS")

GSEA.vacc.other <- GSEA.vacc.other[order(GSEA.vacc.other$SUBGROUP), ]
GSEA.vacc.other$number <- c(1:nrow(GSEA.vacc.other))

GSEA.vacc.other %>% 
  mutate(NAME = fct_reorder(NAME, -number)) %>%
  ggplot( aes(x = sample, y = NAME, label = NAME)) +
  geom_point(aes(size = ABS.NES, colour = Direction, alpha=Significant)) + 
  scale_alpha_discrete(range=c(0.15, 0.6)) + 
  scale_radius(trans="exp", range = c(1,15)) + 
  scale_colour_manual(values=fill) + 
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x=element_blank(), axis.title.y=element_blank()) + 
  ggtitle(label="GSEA on logFC, other modules") +
  theme(panel.background = element_blank())
ggsave("Output/GSEA.on.logFC/GSEA.summary/Plots/Other.modules.pdf", width = 15, height = 12)
