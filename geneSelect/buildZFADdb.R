#####################################################
# ~ ZFAD: build gene database for selection of genes ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# data mainly from Schwartzentruber, 2020
# Nature Genetics doi: 10.1038/s41588-020-00776-w


# packages ----------------------------------------------------------------

library(here)
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)

library(FramebyFrame)


# small functions ---------------------------------------------------------

# find the mode
# from https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# DATA FROM SCHWARTZENTRUBER 2020 -----------------------------------------
# Nature Genetics doi: 10.1038/s41588-020-00776-w

# Table 1: meta-analysis GWAS ---------------------------------------------
zdb <- read_excel(here('geneSelect', 'schwartzentruber2020_suppl.xlsx'), sheet='1-AD loci')

colnames(zdb) <- c('chromosome', 'leadSNP', 'leadSNP_pos', 'leadSNP_freq', 'leadSNP_beta', 'leadSNP_pval',
                   'indSNPs', 'indSNPs_pos', 'indSNPs_n', 'indSNPs_freq', 'indSNPs_beta', 'indSNPs_pval',
                   'locus_start', 'locus_end', 'locus_name', 'nearestgene', 'genes_within100kb', 'genes_within500kb')

# is locus name always same as nearest gene?
zdb[!zdb$locus_name==zdb$nearestgene, c('locus_name', 'nearestgene')]
# 8 cases where it is not the same
# HLA // HLA-DRB1 >> ?
# TREM2 // UNC5CL >> TREM2 almost certainly causal
# IKZF1 // SPATA48 >> ?
# PTK2B-CLU // CLU >> might be both together
# MS4A4A // MS4A6E >> ?
# VKORC1 // KAT8 >> ?
# APP-ADAMTS1 // ADAMTS1 >> APP almost certainly causal


# Table 5: colocalisation -------------------------------------------------
# useful would be:
# one column: which gene is top coloc at this GWAS
# one column: is it the same gene as locus name (i.e. nearest gene is likely causal)
clc <- read_excel(here('geneSelect', 'schwartzentruber2020_suppl.xlsx'), sheet='5-QTL coloc summary') # clc for coloc
colnames(clc) <- c('leadSNP_pos','locus_name', 'signal', 'leadSNP','dataset', 'maxprob','details')

# are locus all the same?
clcloci <- sort(unique(clc$locus_name))
zdbloci <- sort(unique(zdb$locus_name))

# are all coloc loci in zdb loci?
# coloc loci not in zdb loci:
clcloci[!clcloci %in% zdbloci]
# MIR142
# locus name is TSPOAP1 in GWAS/zdb, can see lead SNP is the same
# indeed, TSPOAP1 not in clcloci:
which(clcloci=='TSPOAP1')
# replace MIR142 by TSPOAP1 in coloc
clc[which(clc$locus_name=='MIR142'), 'locus_name'] <- 'TSPOAP1'

# reset clcloci
clcloci <- sort(unique(clc$locus_name))

# try again, are all coloc loci in zdb loci?
unique(clcloci %in% zdbloci) # yes

# in the other way around;
# are all zdb loci in coloc loci?
# zdb loci not in coloc loci:
zdbloci[!zdbloci %in% clcloci] # APOE
# comment from README: "The APOE region was excluded since GWAS p values were too low for reliable coloc."

# from text: taking max coloc score across all tissues works better than composite weighted score based on 'relevant' tissues
# (composite weighted score is like positive colocalisation in microglia has a bigger weight)
# 'works better' because they retrieve a few more genes of a 'truth set' they know from other sources is causal

# will follow this and simply take maximum coloc probability for each locus_name

# details column is list of coloc genes for that dataset, first one is always the highest probability, which is copied in maxprob column
# so need to split that column to keep just first gene name
# format is GENE, probability/...
# so everything before first comma should work
clc$colocgene <- unlist(lapply(strsplit(clc$details, split=','),
                               function(l) {l[1]}))

# preallocate new column in zdb for best coloc gene for each locus
zdb$bestcoloc <- NA
# and probability of that coloc
# because max coloc will always return something but does not say if anything convincing
zdb$bestcoloc_prob <- NA

# for each locus name, get gene with maximum probability
for (l in 1:length(clcloci)) {
  
  # we are looking at locus loc
  loc <- clcloci[l]
  
  # get all the coloc for that gene
  colu <- clc[which(clc$locus_name==loc),] # coloc for that locus
  
  # get the maximum gene and add it to zdb
  zdb[which(zdb$locus_name==loc), 'bestcoloc'] <- 
    as.character(colu[which.max(colu$maxprob), 'colocgene'])
  
  # also get the probability of that coloc and record it in zdb
  zdb[which(zdb$locus_name==loc), 'bestcoloc_prob'] <- 
    as.character(colu[which.max(colu$maxprob), 'maxprob'])
  
}


# -------------------------------------------------------------------------

# alternative, majority vote from datasets where coloc > 0.8
# (highlighted in supplementary table by authors)

# similar approach as above
# then keep only > 0.8 and take statistical mode (i.e. majority vote)
# idea: majority vote from datasets where coloc > 0.8

# preallocate new column in zdb
zdb$majoritycoloc <- NA

# for each locus name, get gene with maximum probability
for (l in 1:length(clcloci)) {
  
  # we are looking at locus loc
  loc <- clcloci[l]
  
  # get all the coloc for that gene
  colu <- clc[which(clc$locus_name==loc),] # coloc for that locus
  
  # keep only the coloc > 0.8
  colu <- subset(colu, maxprob > 0.8)
  
  # find the Mode (majority vote) for that locus and it to zdb
  zdb[which(zdb$locus_name==loc), 'majoritycoloc'] <-
    Mode(colu$colocgene)
  
}


# Table 8: SNP fine-mapping -----------------------------------------------
# from text: missense variants are 19x more likely to be causal; but only 1% of GWAS SNPs
# would be useful here: for each locus, take the variant with highest probability of being causal
  # then highlight if this variant is missense_variant
# may help filling the gaps when coloc is not convincing,
  # i.e. if coloc is not convincing it may simply be because there is a 'simpler explanation',
  # as in, a missense variant is not necessarily predicted to change expression of the gene

# import
fim <- read_excel(here('geneSelect', 'schwartzentruber2020_suppl.xlsx'), sheet='8-SNP Fine-mapping') # fim for fine-mapping

# all we need is the most likely causal variant + probability causal + affected gene (according to VEP, variant effect predictor)

# can do similar to above

# preallocate some new columns in zdb
zdb$topfinemap_snp <- NA
zdb$topfinemap_prob <- NA
zdb$topfinemap_consequence <- NA
zdb$topfinemap_gene <- NA
zdb$topfinemap_biotype <- NA

# are all loci present in fine-mapping sheet?
fimloci <- sort(unique(fim$`locus name`))
zdbloci <- sort(unique(zdb$locus_name))

# are all fine-map loci in zdb loci?
# fine-map loci not in zdb loci:
fimloci[!fimloci %in% zdbloci]
# yes

# opposite: are all zdb loci in fine-map loci?
# zdb loci not in fine-map loci:
zdbloci[!zdbloci %in% fimloci] # APOE
# must be same reason as above

# for each locus name, get SNP with maximum probability
for (l in 1:length(fimloci)) {
  
  # we are looking at locus loc
  loc <- fimloci[l]
  
  # get all the fine-mapping SNPs for that gene
  snps <- fim[which(fim$`locus name`==loc),] # fine-mapping SNPs for that locus
  
  # which one is the best according to mean probability?
  toprow <- which.max(snps$`mean prob`) # row of top SNP
  
  # fill in zdb columns with info for that SNP
  zdb[which(zdb$locus_name==loc), 'topfinemap_snp'] <- snps[toprow, 'snp']
  zdb[which(zdb$locus_name==loc), 'topfinemap_prob'] <- snps[toprow, 'mean prob']
  zdb[which(zdb$locus_name==loc), 'topfinemap_consequence'] <- snps[toprow, 'Consequence']
  zdb[which(zdb$locus_name==loc), 'topfinemap_gene'] <- snps[toprow, 'SYMBOL']
  zdb[which(zdb$locus_name==loc), 'topfinemap_biotype'] <- snps[toprow, 'BIOTYPE']
  
}



# Table 10: network ranking ------------------------------------------------
# from text: degree to which the gene is supported by its interaction with top AD candidate genes across all other loci
# does not give info per locus, but for pretty much every gene (18k rows)
# network score is copied for genes to be prioritised in Table 13
# so will simply take from there


# Table 13: summary prioritisation ----------------------------------------
# most likely causal gene from five lines of evidence

# import
prio <- read_excel(here('geneSelect', 'schwartzentruber2020_suppl.xlsx'), sheet='13-gene evidence rankings')

# as before, preallocate a few columns in zdb

# I think column model prob is essentially same as totalscore, just brought down to 0--1
# but will keep both
zdb$modelprob_gene <- NA # top gene by model probability
zdb$modelprob <- NA # model probability for that gene

zdb$totalscore_gene <- NA # top gene by total score (likely same as by model probability)
zdb$totalscore_geneID <- NA # same, but Ensembl ID
zdb$totalscore <- NA # total score for that gene

# total score is sum coding score + coloc score + distance score + expression score (all taken below, except coloc)

zdb$distscore_gene <- NA # top gene by distance score, i.e. how far is potentially causal gene from lead GWAS SNP
zdb$distscore <- NA # distance score for that gene

zdb$codingscore_gene <- NA # top gene by coding score; most of information here was probably extracted when looked at fine-mapping genes
zdb$codingscore <- NA # coding score for that gene

# skipping coloc score column, we extracted that information above 

zdb$networkscore_gene <- NA # top gene by network score, see notes above
zdb$networkscore <- NA # network score for that gene

zdb$exprscore_gene <- NA # expression score for specificity of expression in microglia or brain, whichever is max
zdb$exprscore <- NA # expression score for that gene

# similar approach as other sections

# first, are all loci there?
prioloci <- sort(unique(prio$locus))
identical(prioloci, zdbloci) # yes

# for each locus name, find the best gene by total score etc.
for (l in 1:length(prioloci)) {
  
  # we are looking at locus loc
  loc <- prioloci[l]
  
  # get the information for all the possibly causal genes of that locus
  gens <- prio[which(prio$locus==loc),]
  
  # by modelprob
  # which gene is best?
  zdb[which(zdb$locus_name==loc), 'modelprob_gene'] <- gens[which.max(gens$`model prob`), 'symbol']
  # model prob of that gene?
  zdb[which(zdb$locus_name==loc), 'modelprob'] <- gens[which.max(gens$`model prob`), 'model prob']
  
  # by totalscore
  # which gene is best?
  zdb[which(zdb$locus_name==loc), 'totalscore_gene'] <- gens[which.max(gens$`total score`), 'symbol']
  # Ensembl ID of that gene
  zdb[which(zdb$locus_name==loc), 'totalscore_geneID'] <- gens[which.max(gens$`total score`), 'geneID']
  # total score of that gene?
  zdb[which(zdb$locus_name==loc), 'totalscore'] <- gens[which.max(gens$`total score`), 'total score']
  
  # by distscore
  # which gene is best?
  zdb[which(zdb$locus_name==loc), 'distscore_gene'] <- gens[which.max(gens$`gene dist score`), 'symbol']
  # distance score of that gene?
  zdb[which(zdb$locus_name==loc), 'distscore'] <- gens[which.max(gens$`gene dist score`), 'gene dist score']
  
  # by codingscore
  # which gene is best?
  zdb[which(zdb$locus_name==loc), 'codingscore_gene'] <- gens[which.max(gens$`coding score`), 'symbol']
  # coding score of that gene?
  zdb[which(zdb$locus_name==loc), 'codingscore'] <- gens[which.max(gens$`coding score`), 'coding score']
  
  # by network score
  # which gene is best?
  zdb[which(zdb$locus_name==loc), 'networkscore_gene'] <- gens[which.max(gens$`network score`), 'symbol']
  # coding score of that gene?
  zdb[which(zdb$locus_name==loc), 'networkscore'] <- gens[which.max(gens$`network score`), 'network score']
  
  # by expression score
  # which gene is best?
  zdb[which(zdb$locus_name==loc), 'exprscore_gene'] <- gens[which.max(gens$`expr score`), 'symbol']
  # expression score of that gene?
  zdb[which(zdb$locus_name==loc), 'exprscore'] <- gens[which.max(gens$`expr score`), 'expr score']
  
}


### done with extracting info from SCHWARTZENTRUBER 2020


# some summary info -------------------------------------------------------

# when is the nearest gene NOT the most likely causal gene?
zdb[!zdb$nearestgene == zdb$totalscore_gene, c('nearestgene', 'totalscore_gene')]
# 9 cases, out of 37
# includes
# HLA-DRB1 vs HLA-DRA
# so really 8 cases

# locus name and nearest gene can be a little bit different, so:
# when is the locus name NOT the most likely causal gene?
zdb[!zdb$locus_name == zdb$totalscore_gene, c('locus_name', 'totalscore_gene')]
# 10 cases, out of 37
# includes
# HLA vs HLA-DRA
# PTK2B-CLU vs PTK2B
# so really 8 cases, out of 37


# partial summary causal gene ---------------------------------------------

# will make a new column causal gene (most likely causal gene) & its Ensembl ID
# will put best by total score for now, then correct some manually
zdb$causalgene <- zdb$totalscore_gene
zdb$causalgeneID <- zdb$totalscore_geneID

# what we can correct already:
# we know causal gene at TREM2 is TREM2 from rare variants, not TREM1
zdb[which(zdb$locus_name=='TREM2'), 'causalgene'] <- 'TREM2'
zdb[which(zdb$locus_name=='TREM2'), 'causalgeneID'] <- 'ENSG00000095970'

# at PTK2B-CLU locus, current causal gene is PTK2B
# but there is evidence that they may both be interesting. May even interact
# will do two unique rows, one for CLU one for PTK2B. Most entries will be the same
zdb <- rbind(zdb, zdb[which(zdb$locus_name=='PTK2B-CLU'),]) # copy the row
zdb[nrow(zdb), 'causalgene'] <- 'CLU' # change causal gene to CLU
zdb[nrow(zdb), 'causalgeneID'] <- 'ENSG00000120885'

# at APP-ADAMTS1 locus, we know causal gene is APP from trisomies, amyloid beta etc.
zdb[which(zdb$locus_name=='APP-ADAMTS1'), 'causalgene'] <- 'APP'
zdb[which(zdb$locus_name=='APP-ADAMTS1'), 'causalgeneID'] <- 'ENSG00000142192'


# deal with familial AD genes ---------------------------------------------
# will be useful to add them now so they get integrated in the subsequent analyses
# APP is already there (APP-ADAMTS1 locus), but not PSEN1 and PSEN2

# add PSEN1
zdb <- rbind(zdb, NA)
zdb[nrow(zdb), 'chromosome'] <- 14
zdb[nrow(zdb), 'locus_name'] <- 'PSEN1'
zdb[nrow(zdb), 'causalgene'] <- 'PSEN1'
zdb[nrow(zdb), 'causalgeneID'] <- 'ENSG00000080815'

# add PSEN2
zdb <- rbind(zdb, NA)
zdb[nrow(zdb), 'chromosome'] <- 1
zdb[nrow(zdb), 'locus_name'] <- 'PSEN2'
zdb[nrow(zdb), 'causalgene'] <- 'PSEN2'
zdb[nrow(zdb), 'causalgeneID'] <- 'ENSG00000143801'


# DATA FROM RAJ 2018 ------------------------------------------------------
# in Nature Genetics: doi.org/10.1038/s41588-018-0238-1

# ! I think uses hg19 throughout

# Table 2: spliced introns w/ neuropathologies ----------------------------
# list of significantly differentially spliced introns associated with neuropathologies
# here: only interested if it corroborates GWAS gene, as gives more evidence of causation
# e.g. most significant differentially spliced intron is in PFKP, but it is not a GWAS locus, so could be consequence of disease
# differentially spliced intron in gene can be associated with neuritic plaques / amyloid / tangles

sin <- read_excel(here('geneSelect', 'raj2018_suppl', 'raj2018_suppltable2.xlsx')) # sin for spliced intro vs neuropathologies

subset(sin,
       gene_id %in% zdb$causalgene) # only one is APP, so not incredibly useful

# can add that info to zdb anyways

# inspliPatho = intron spliced w/ pathology
inspliPatho <- subset(sin, gene_id %in% zdb$causalgene, c('gene_id', 'Z-score', 'P-value', 'Trait'))

colnames(inspliPatho) <- c('splintron_neurop_geneID', 'splintron_neurop_zscore', 'splintron_neurop_pval', 'splintron_neurop_trait')

# add to zdb
zdb <- merge(zdb, inspliPatho, by.x='causalgene', by.y='splintron_neurop_geneID', all.x=TRUE)


# Table 3: spliced introns w/ diagnosis -----------------------------------

sid <- read_excel(here('geneSelect', 'raj2018_suppl', 'raj2018_suppltable3.xlsx')) # sid for spliced intron with diagnosis

# ! this list includes all (or most) genes (27681 rows), but many have adjusted p-value = 1
# will do as in text and use Bonferonni pval < 0.05
# inspliDiagno = intron spliced w/ diagnosis
inspliDiagno <- subset(sid,
                       gene_id %in% zdb$causalgene &
                         `P-value_BonferroniAdjusted` < 0.05,
                       c('gene_id', 'P-value_BonferroniAdjusted'))
# CLU, APP, PICALM

colnames(inspliDiagno) <- c('splintron_diagno_geneID', 'splintron_diagno_pvalBonferonni')

# add to zdb
zdb <- merge(zdb, inspliDiagno, by.x='causalgene', by.y='splintron_diagno_geneID', all.x=TRUE)


# Table 8: replications of spliced introns w/ diagnosis -------------------

# will skip here to keep it simple
# but important to note: of CLU, APP, PICALM, only APP is replicated
# dataset is a bit smaller though, which may be the reason why it does not come up


# Table 9: spliced introns w/ overexpression Tau --------------------------
# list of differentially spliced introns from overexpressing Tau in iPSC-derived neurons
sit <- read_excel(here('geneSelect', 'raj2018_suppl', 'raj2018_suppltable9.xlsx')) # sit for spliced intron with Tau overexpression

# ! most/all genes included (21,650 rows) so need to use p-val
# inspliTau = intron spliced w/ Tau overexpression in iPSC
inspliTau <- subset(sit,
                    gene %in% zdb$causalgene &
                      `P-value_Benjamini-Hochberg` < 0.05,
                    c('gene', 'P-value_Benjamini-Hochberg'))
# APP, PICALM

colnames(inspliTau) <- c('splintron_tau_geneID', 'splintron_tau_pval_BenjaminiHochberg')

# add to zdb
zdb <- merge(zdb, inspliTau, by.x='causalgene', by.y='splintron_tau_geneID', all.x=TRUE)


# Table 11: TWAS ----------------------------------------------------------
twas <- read_excel(here('geneSelect', 'raj2018_suppl', 'raj2018_suppltable11.xlsx'))

# TWAS: splicing
twas_spli <- twas[2:17,]

twas_spli <- subset(twas_spli,
                    ID %in% zdb$causalgene,
                    c('ID', 'TWAS.Z', 'TWAS.P'))
colnames(twas_spli) <- c('TWASspli_geneID', 'TWASspli_zscore','TWASspli_pval')
# CLU, PTK2B, PICALM

# add to zdb
zdb <- merge(zdb, twas_spli, by.x='causalgene', by.y='TWASspli_geneID', all.x=TRUE)

###

# TWAS: gene expression
twas_exp <- twas[20:24,]

twas_exp <- subset(twas_exp,
                   ID %in% zdb$causalgene,
                   c('ID', 'TWAS.Z', 'TWAS.P'))
colnames(twas_exp) <- c('TWASexp_geneID', 'TWASexp_zscore','TWASexp_pval')
# CR1, SPI1

# add to zdb
zdb <- merge(zdb, twas_exp, by.x='causalgene', by.y='TWASexp_geneID', all.x=TRUE)


# Table 12: second TWAS ---------------------------------------------------
# I think it is splicing, not expression level

twas2 <- read_excel(here('geneSelect', 'raj2018_suppl', 'raj2018_suppltable12.xlsx'))

twas2 <- subset(twas2,
                ID %in% zdb$causalgene,
                c('ID', 'TWAS.Z', 'TWAS.P'))
colnames(twas2) <- c('TWASspli2_geneID', 'TWASspli2_zscore','TWASspli2_pval')
# ABCA7

# add to zdb
zdb <- merge(zdb, twas2, by.x='causalgene', by.y='TWASspli2_geneID', all.x=TRUE)


# summary from Raj 2018 ---------------------------------------------------

# what did we learn?
# APP is supported by intron splicing; APP was already set as causal gene at APP-ADAMTS1 locus
# both CLU and PTK2B are supported by intron splicing; both already set as causal genes at the CLU-PTK2B locus
# PICALM supported by intron splicing; PICALM was already set as causal gene at PICALM locus
# SPI supported by (expression) TWAS; SPI1 was already set as causal gene at SPI1 locus

# so does not change any causal gene but extra support helps with prioritisation


# ZEBRAFISH ORTHOLOGUES ---------------------------------------------------
# of most likely causal gene at each locus

# how I proceeded;

# How to get human/zebrafish Ensembl orthologs

# Essentially following http://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
  
# http://www.ensembl.org/biomart/martview
# 'Ensembl Genes 109'
# Human genes (hg38)
# Filters
# Multi Species Comparisons > tick Homologues filters > Orthologous Zebrafish Genes
# Attributes > Homologues > Orthologues [U-Z] > choose which attributes wanted under Zebrafish Orthologues (tick all)
# Results (top)
# Export all results to File CSV
# saves in .txt but it is a CSV, change name to
# zebrafishOrthologues.csv

# last updated 27/04/2023

# import
ortho <- read.csv(here('geneSelect', 'zebrafishOrthologues.csv'))

# some format conversion issues
ortho$Gene.stable.ID <- as.character(ortho$Gene.stable.ID)
ortho$Gene.stable.ID.version <- as.character(ortho$Gene.stable.ID.version)
ortho$Transcript.stable.ID <- as.character(ortho$Transcript.stable.ID)
ortho$Transcript.stable.ID.version <- as.character(ortho$Transcript.stable.ID.version)
ortho$Zebrafish.gene.stable.ID <- as.character(ortho$Zebrafish.gene.stable.ID)
ortho$Zebrafish.gene.name <- as.character(ortho$Zebrafish.gene.name)
ortho$Zebrafish.protein.or.transcript.stable.ID <- as.character(ortho$Zebrafish.protein.or.transcript.stable.ID)
ortho$Zebrafish.chromosome.scaffold.name <- as.character(ortho$Zebrafish.chromosome.scaffold.name)
ortho$Zebrafish.chromosome.scaffold.start..bp. <- as.numeric(ortho$Zebrafish.chromosome.scaffold.start..bp.)
ortho$Zebrafish.chromosome.scaffold.end..bp. <- as.numeric(ortho$Zebrafish.chromosome.scaffold.end..bp.)
ortho$Query.protein.or.transcript.ID <- as.character(ortho$Query.protein.or.transcript.ID)
ortho$Last.common.ancestor.with.Zebrafish <- as.character(ortho$Last.common.ancestor.with.Zebrafish)
ortho$Zebrafish.homology.type <- as.character(ortho$Zebrafish.homology.type)
ortho$X.id..target.Zebrafish.gene.identical.to.query.gene <- as.numeric(ortho$X.id..target.Zebrafish.gene.identical.to.query.gene)
ortho$X.id..query.gene.identical.to.target.Zebrafish.gene <- as.numeric(ortho$X.id..query.gene.identical.to.target.Zebrafish.gene)
ortho$Zebrafish.Gene.order.conservation.score <- as.numeric(ortho$Zebrafish.Gene.order.conservation.score)
ortho$Zebrafish.Whole.genome.alignment.coverage <- as.numeric(ortho$Zebrafish.Whole.genome.alignment.coverage)
ortho$Zebrafish.orthology.confidence..0.low..1.high. <- as.logical(ortho$Zebrafish.orthology.confidence..0.low..1.high.)

# about columns
# geneID
# geneID + version
# transcript ID
# transcript ID + version
# zebrafish gene ID
# zebrafish gene name
# zebrafish protein or transcript ID
# zebrafish chromosome
# zebrafish start position (bp)
# zebrafish end position (bp)
# query protein or transcript ID
# last common ancestor with zebrafish
# homology type
  # one2one = one human gene vs one zebrafish gene
  # one2many = one human gene vs multiple zebrafish genes (i.e. gene underwent duplication in zebrafish lineage, after common ancestor)
  # many2many = gene underwent duplication in both lineage, after common ancestor
# % identity target gene (zebrafish) vs query gene (human)
# % identity query gene (human) vs target gene (zebrafish)
  # not sure why would be different; and they are very similar
  # I guess if lengths of gene are different; need to decide which is reference?
# zebrafish gene-order conservation score
  # see https://m.ensembl.org/info/genome/compara/Ortholog_qc_manual.html
  # basically: take the two genes before and the two after in each species
  # do they match between species (i.e. orthologues, and in same order)
  # each is worth 25%
# zebrafish whole-genome alignment coverage
  # see https://m.ensembl.org/info/genome/compara/Ortholog_qc_manual.html
  # basically how well the two genes (genomic window) align together
  # exon aligning is considered more important than intron aligning (but bonus points if introns align)
# zebrafish orthology confidence (0 or 1 converted to FALSE/TRUE)
  # 1 = tagged as high confidence
  # see https://m.ensembl.org/info/genome/compara/Ortholog_qc_manual.html
  # just follows some thresholds for % identity / gene-order conservation score / whole-genome alignment coverage
  # threshold depends on common ancestor: if more recent, thresholds to meet are higher

# will need to match with Ensembl ID
# i.e. column causalgeneID from zdb
zdb$causalgeneID %in% ortho$Gene.stable.ID # not all are in Ensembl orthology database
# which ones are not in Ensembl orthology database?
zdb[!zdb$causalgeneID %in% ortho$Gene.stable.ID, c('causalgene', 'causalgeneID')]

# for those;
# checked all Ensembl ID; all looks correct
# check if any gene with that name in ZFIN.org

# Note;
# ADAM10 had no ortholog recorded in 2021, now (2022) it has some entries
ortho[which(ortho$Gene.stable.ID=='ENSG00000137845'),]

# CD33 had no ortholog recorded in 2021, now (2022) it has some entries
ortho[which(ortho$Gene.stable.ID=='ENSG00000105383'),]
# update 2023: lost its orthologue(s)?

# still without any entry:
# ACE
# CR1
# EPHA1
# HLA-DRA
# HS3TS1
# PILRA
# SCIMP
# SPPL2A
# TREM2

# update 2023: CD33 lost its orthologue(s)?
# no difference for other genes

# all Ensembl ID seem correct/up to date
# so means no orthologue found for these genes in Ensembl database

# checking zebrafish genes with same names in Ensembl & ZFIN
# ace exists; ENSDARG00000079166
# cr1 does not
# epha1 does not
# hla-dra (or hla alone) does not
# hs3ts1 does not
# pilra does not
# scimp does not
# sppl2a does not
# trem2 does not

# >> need to investigate ACE and potentially add it manually

# often multiple entries for orthologues for one human gene
# so make a list of dataframes, each dataframe = orthologues of one human gene
zfOrths <- vector(mode='list', length=nrow(zdb))
names(zfOrths) <- c(zdb$causalgene)

# fill in the list
for (h in 1:length(zdb$causalgene)) {
  hid <- zdb$causalgeneID[h] # looking for orthologues of that human gene ID
  
  # make a small temporary dataframe with the results
  oth <- ortho %>%
    subset(Gene.stable.ID==hid) %>%
    add_column(humangene=names(zfOrths)[[h]], .before='Gene.stable.ID') # add human gene name as first column
  
  # add that dataframe in the list
  zfOrths[[h]] <- oth
  
}

# now summarise the list
# keep just one row per zebrafish gene

for (h in 1:length(zfOrths)) {
  # each element of the list is one human gene (hence h)
  # what are the unique zebrafish orthologues for that human gene?
  zids <- unique(zfOrths[[h]]$Zebrafish.gene.stable.ID) # zebrafish gene IDs
  
  # only keep one entry for each zebrafish gene ID
  # so preallocate a small temporary dataframe, one row per zebrafish gene ID
  zos <- as.data.frame(matrix(nrow=length(zids), ncol=ncol(zfOrths[[h]]))) # zebrafish orthologues
  colnames(zos) <- colnames(zfOrths[[h]])
    
  # now fill in one row for each zebrafish gene ID
  for (z in 1:length(zids)) {
    zos[z,] <- zfOrths[[h]] [which(zfOrths[[h]]$Zebrafish.gene.stable.ID==zids[z])[1] , ] # which...[1] finds the first entry (row index) for that zid
  }
  
  # add the human gene name in zos
  # if some orthologue, human gene name will be transferred above
  # but if no orthologue, it is now NA
  # so puts back the human gene name either way so we always have it
  zos$humangene <- names(zfOrths)[h]
  
  # zos for that human gene is ready, now we replace corresponding element in list zfOrths
  zfOrths[[h]] <- zos
}


###

# clean up in a nicer dataframe we can export
# will do a separate one than causal gene to not add too many columns when more than one orthologue
# one row per zebrafish orthologue
# simply rbind the list
zfos <- rbindlist(zfOrths)
# it puts human genes

# clean up the columns
# some columns we do not need
zfos$Gene.stable.ID.version <- NULL
zfos$Transcript.stable.ID <- NULL
zfos$Transcript.stable.ID.version <- NULL

# update column names
colnames(zfos) <- c('humangene', 'hgeneID', 'zgeneID', 'zgene', 'zproteinID', 'zchr', 'zendpos', 'zstartpos', 'hproteinID',
                    'lastcommonancestor', 'homologytype', 'geneidentity_z2h', 'geneidentity_h2z', 'geneOrderScore', 'alignmentScore', 'highConfidence')



# zdb: order columns ------------------------------------------------------

# delete a bunch of not-very-useful columns
zdb$indSNPs <- NULL
zdb$indSNPs_pos <- NULL
zdb$indSNPs_n <- NULL
zdb$indSNPs_freq <- NULL
zdb$indSNPs_beta <- NULL
zdb$indSNPs_pval <- NULL
zdb$totalscore_geneID <- NULL

# rename a few columns
colnames(zdb)[which(colnames(zdb)=='bestcoloc')] <- 'bestcoloc_gene'
colnames(zdb)[which(colnames(zdb)=='majoritycoloc')] <- 'majoritycoloc_gene'

# reorder columns
newcolorder <- c('locus_name', 'causalgene', 'causalgeneID', 'chromosome', 'locus_start', 'locus_end',
                 'leadSNP', 'leadSNP_pos', 'leadSNP_freq', 'leadSNP_beta', 'leadSNP_pval', 'nearestgene',
                 'genes_within100kb', 'genes_within500kb', 'bestcoloc_gene', 'bestcoloc_prob', 'majoritycoloc_gene',
                 'topfinemap_snp', 'topfinemap_prob', 'topfinemap_consequence', 'topfinemap_gene', 'topfinemap_biotype',
                 'distscore_gene', 'distscore', 'codingscore_gene', 'codingscore', 'networkscore_gene', 'networkscore',
                 'exprscore_gene', 'exprscore', 'modelprob_gene', 'modelprob', 'totalscore_gene', 'totalscore',
                 'splintron_neurop_zscore', 'splintron_neurop_pval', 'splintron_neurop_trait', 'splintron_diagno_pvalBonferonni',
                 'splintron_tau_pval_BenjaminiHochberg', 'TWASspli_zscore', 'TWASspli_pval', 'TWASexp_zscore', 'TWASexp_pval', 'TWASspli2_zscore', 'TWASspli2_pval')

# check all columns are there
unique(colnames(zdb) %in% newcolorder)

# now put zdb in new order
zdb <- zdb[match(newcolorder, colnames(zdb))]



# zfos: order columns -----------------------------------------------------

zfos <- as.data.frame(zfos)

colnames(zfos)[which(colnames(zfos)=='zgene')] <- 'zebrafishgene'

newcolorder <- c('humangene', 'hgeneID', 'hproteinID', 'zebrafishgene', 'zgeneID', 'zproteinID', 
                 'zchr', 'zstartpos', 'zendpos', 'lastcommonancestor', 'homologytype',
                 'geneidentity_z2h', 'geneidentity_h2z', 'geneOrderScore', 'alignmentScore', 'highConfidence')

# check all columns are there
unique(colnames(zfos) %in% newcolorder)

# now put zfos in new order
zfos <- zfos[match(newcolorder, colnames(zfos))]



# some more edits on orthologues ------------------------------------------

# column homologytype is written e.g. ortholog_one2one
# just keep one2one
zfos$homologytype <- strNthSplit(zfos$homologytype, '_', 2)


### some odd results

# zebrafish ABCA7 is written like a human gene; it seems to be written like this in Ensembl
# https://www.ensembl.org/Danio_rerio/Gene/Summary?db=core;g=ENSDARG00000074221;r=11:5588122-5659681;t=ENSDART00000113281
# ZFIN https://zfin.org/ZDB-GENE-050517-6 writes it as a zebrafish gene but links it to same Ensembl page, so will leave as it is

# zebrafish aph1b, written like an ohnologue but I do not see any aph1a zebrafish gene
# consistent with the fact that it is orthologue of APH1B, not APH1

# surprising that orthologue of SPI1 is spi1b, and no mention of spi1a
# why not spi1a?
# ZFIN says spi1a is orthologue of SPI1
# https://zfin.org/ZDB-GENE-060825-351#orthology
# comparing amino acid sequence:

# Zebrafish protein spi1a
# MLHYRMESCVISPLSEEIIPYEHEARPIYDFYPYLSTDPETHPESGCEYSSVYGHHSEFE
# PPPGSHFTELHTSTYRYGDMEAFHPGVDAAMGTILPVVPPQYTYISHPLYQRSPVPHCST
# DEEEPGGRSPPFEVSEGEEDHDGHPSTSSTLSGNKRKVRLYQFLLDLLQDGDMRDCIWWV
# DRERGVFQFSSKHKETLASRWGQQKGNRKRMTYQKMARALRNYGKTGEVKKVKKKLTYQF
# SGDVLRRVTMERRQYHH

# Human protein SPI1
# MLQACKMEGFPLVPPPSEDLVPYDTDLYQRQTHEYYPYLSSDGESHSDHYWDFHPHHVHS
# EFESFAENNFTELQSVQPPQLQQLYRHMELEQMHVLDTPMVPPHPSLGHQVSYLPRMCLQ
# YPSLSPAQPSSDEEEGERQSPPLEVSDGEADGLEPGPGLLPGETGSKKKIRLYQFLLDLL
# RSGDMKDSIWWVDKDKGTFQFSSKHKEALAHRWGIQKGNRKKMTYQKMARALRNYGKTGE
# VKKVKKKLTYQFSGEVLGRGGLAERRHPPH

# >> 50% identity

# Zebrafish protein spi1b
# MLHPYRMEGYIIPPKEEGRERVTWTGWMSQTPSVQKDYWAVLTKDQQTEEMFETEIYRPP
# MEYQYIIDDSQNDHSWDYNTHHIHPVDFENLPESHFTELQSVQSLHAASVHRFPDVESSH
# FMDPGLGSHHIPLATPQMTYLPRTSVCYPHNVQPSPLQRSSDEEDPSSRSPPLEVSDEEC
# MRDHISSTTGGEHGNKKKIRLYQFLLDLLRNGDMKDSIWWVDREKGTFQFSSKHKEVLAN
# RWGIQKGNRKKMTYQKMARALRNYGKTGEVKKIKKKLTYQFSGEVLGKSHTDRKHYM

# vs SPI1 human >> 61% identity

# looks OK I think? 50% identity is higher than many
# so I do not know why spi1a not reported as orthologue by Ensembl
# will add it manually

zfos <- rbind(zfos, data.frame(humangene='SPI1', hgeneID='ENSG00000066336', hproteinID='ENSP00000367799',
                               zebrafishgene='spi1a', zgeneID='ENSDARG00000067797', zproteinID='ENSDARP00000098095',
                               zchr=25, zstartpos=35553542, zendpos=-35572245, lastcommonancestor=NA, homologytype=NA,
                               geneidentity_z2h=NA, geneidentity_h2z=NA, geneOrderScore=NA, alignmentScore=NA, highConfidence=NA))


# write zdb & zfos as two separate sheets of same file --------------------

zfwb <- createWorkbook() # ZF workbook

modifyBaseFont(zfwb, fontSize=11, fontColour="black", fontName="Courier New")

headerstyle <- createStyle(fontName='Courier New', fontSize=11, fontColour='#1F497D',
                           textDecoration='bold', border='bottom', borderColour='#1F497D', borderStyle='medium')

addWorksheet(zfwb, 'gwas_causalgene')
addWorksheet(zfwb, 'zebrafish_orthologues')

writeData(zfwb, sheet='gwas_causalgene', zdb)
writeData(zfwb, sheet='zebrafish_orthologues', zfos)

addStyle(zfwb, 'gwas_causalgene', headerstyle, rows=1, cols=1:ncol(zdb))
addStyle(zfwb, 'zebrafish_orthologues', headerstyle, rows=1, cols=1:ncol(zfos))

# highlight in zdb all which do not have a zebrafish orthologue
# which gene does not have a zebrafish orthologue?
noorth <- zfos[which(is.na(zfos$zebrafishgene)), 'humangene'] # human genes without an orthologue

# create a style for those, light red highlight & greyed out font
noOrthstyle <- createStyle(fontName='Courier New', fontSize=11, fontColour='#404040', fgFill='#F2DCDB')

# apply this style to the rows of causal genes without orthologue
addStyle(zfwb, 'gwas_causalgene', noOrthstyle, rows=match(noorth, zdb$causalgene)+1, cols=1:ncol(zdb), gridExpand=TRUE)

# also highlight genes without orthologue in sheet zebrafish_orthologues
addStyle(zfwb, 'zebrafish_orthologues', noOrthstyle, rows=match(noorth, zfos$humangene)+1, cols=1:ncol(zfos), gridExpand=TRUE)

# write the workbook
saveWorkbook(zfwb, file=here('geneSelect', 'ZFAD.xlsx'), overwrite=TRUE)
