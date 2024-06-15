#####################################################
# ~ ZFAD: exploring single-cell RNAseq data from Raj et al., 2020 (Schier lab)
# looking at zebrafish orthologues of Alzheimer's disease genes ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################


# packages ----------------------------------------------------------------

library(here)
library(RCurl)
library(openxlsx)
# BiocManager::install("GEOquery")
library(GEOquery)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(tibble)
library(data.table)
library(ggbeeswarm)
library(gtools)

library(FramebyFrame)



# quick tutorial on 6 somites data ----------------------------------------

# note for @counts and others to work, I had to reinstall a package called Matrix

# sam for sample
sam <- readRDS(here('zfBrain_scRNAseq', 'GSE158142_zf6s_cc_filt.cluster.rds'))
sam <- UpdateSeuratObject(sam) # for the big list seuratAll, this was done in prepData.R

sam[['RNA']]@counts[1:100, 1:10] # I think this is columns = cells / rows = genes / data = number of transcripts
sam[['RNA']]@data[1:100, 1:10] # I think this is columns = cells / rows = genes / data = ?; probably expression level

# we could tell for sure
# e.g. apoea, there are two cells > 0.7
apoea <- sam[['RNA']]@data [which(rownames(sam[['RNA']]@data)=='apoea'),]
apoea[which(apoea > 0.7)] # yes, two cells ~ 0.8, which matches plot

# this way we can put a threshold e.g. minimum 3 cells with expression

# where do we get cluster assignment?
# one of the $res.X, but I do not know which one
# we can find this way
length(unique(sam$res.1))
length(unique(sam$res.1.5))
length(unique(sam$res.2))
length(unique(sam$res.2.5))
length(unique(sam$res.3))
length(unique(sam$res.3.5))
length(unique(sam$res.4))
length(unique(sam$res.4.5))
length(unique(sam$res.5))
# and zfclusters.csv for 12 hpf has 44 clusters
# which does not correspond to any...

# Methods do not say but say as previously described in
# https://www.nature.com/articles/nbt.4103#Sec9
# below Cell type clustering analysis., talks about resolution 2.5

length(sort(as.numeric(unique(sam$res.2.5))))
# gives 38, numbers from 0 to 38
# numbers in zfclusters.csv also give 0 to 45
# hard to understand...

# ! clustering resolution listed here
# https://github.com/brlauuu/zf_brain/blob/master/data/stage_time_clustering_res_mapping.tsv
# so 6 somites is 4.5
# which is 46 clusters
sort(as.numeric(unique(sam$res.4.5)))
# 0 to 45
# OK I understand, some clusters are missing from zfclusters.csv
# cluster3
# cluster6
# so 44 unique clusters in zfclusters.csv + 2 missing = 46 clusters
# I guess unidentified?


# import ------------------------------------------------------------------

# note, we only need all the data to generate expressed/not expressed grid plot
# avoid importing every time as it takes time
seus <- readRDS(here('zfBrain_scRNAseq', 'seuratAll.rds'))

# if only need to update violin plots, can import the 5-dpf data only
# seu5dpf <- readRDS(here('zfBrain_scRNAseq', 'seurat5dpf.rds'))

cls <- read.csv(here('zfBrain_scRNAseq', 'zfclusters.csv'))

# in Seurat data, clusters are as character
# do the same in cls
cls$cluster <- as.character(cls$cluster)


# plot for one gene/one timepoint -----------------------------------------
# e.g. apoeb

seus[[1]]

VlnPlot(seus[[1]], features='apoeb') # normalised counts
VlnPlot(seus[[1]], features='apoeb', slot='counts') # raw counts

FeaturePlot(seus[[1]], features='apoeb', label=TRUE)

# add cluster information


# what is 28
subset(cls,
       stage=='12 hpf' & cluster==28)

subset(cls,
       stage=='12 hpf' & cluster==1)


# plot for one gene/all timepoints ----------------------------------------

# polished violin plot for apoeb

gene <- 'apoeb'

# list to store ggplots
vlns <- vector(length=length(seus), mode='list') # violin plots

# fill it with the violin plots
for (i in 1:length(seus)) {
  vlns[[i]] <- 
    VlnPlot(seus[[i]], features=gene) +
    ylab('expression') +
    ggtitle(names(seus)[[i]]) +
    theme(
      legend.position='none',
      axis.title.x=element_blank()
    )
}

# arrange in a grid
vlnsGrid <- ggarrange(plotlist=vlns, nrow=length(vlns), ncol=1)

# save to a pdf
dir.create(here('zfBrain_scRNAseq', 'grids'))
ggsave(filename=here('zfBrain_scRNAseq', 'grids', paste0(gene, '.pdf')), plot=vlnsGrid, width=400, height=1000, units='mm')



# plot grid there/not there -----------------------------------------------
# grid where each row = one developmental stage, each column = one AD gene

# import list of AD genes (precisely: orthologues of human AD genes)
orth <- read.xlsx(here('geneSelect', 'ZFAD.xlsx'), sheet='zebrafish_orthologues')

# delete the orthologues of MS4A6E, too many of them
orth <- orth[which(orth$humangene!='MS4A6E'),]
# so AD genes we are interested in are
adgenes <- sort(unique(orth$zebrafishgene))

# note, ABCA7 seems to be correct nomenclature in Ensembl
# but previous name is zgc:172302, which seems to be present in data
# so change in adgenes so finds it correctly
adgenes[which(adgenes=='ABCA7')] <- 'zgc:172302'

# other that is not found at any developmental stage is CABZ01076737.1'
# it is orthologue of TSPOAP1
# there seems to be a tspoap1 gene on ZFIN: https://zfin.org/ZDB-GENE-161021-1#summary
# but also does not exist in dataset
# adgenes[which(adgenes=='CABZ01076737.1')] <- 'tspoap1'
# so will leave it as it is


# are genes expressed at the different developmental stage?
thno <- lapply(1:length(names(seus)), function(st) {
  
  sapply(adgenes, function(g) {
    g %in% rownames(seus[[st]])
  })
  
})


### v2, minimum 3 cells (i.e. 4 or more) with expression, makes no difference
# minimum 10 cells makes some difference, fyi
thno <- lapply(1:length(seus), function(st) {
  # loop through developmental stages
  cat('\n \t \t \t \t >>> Developmental stage:', names(seus)[st], '\n')
  
  sapply(adgenes, function(g) {
    # loop through genes
    cat('\t \t \t \t \t >>>', g, '\n')
    
    # is the gene even present?
    # if not, return FALSE, i.e. absent
    if( ! g %in% rownames(seus[[st]][['RNA']]@data) ) {
      cat('\t \t \t \t \t \t >>> not present \n')
      return(FALSE)
    }
    
    # if yes, get expression of every cell
    excells <- seus[[st]][['RNA']]@data [which(rownames(seus[[st]][['RNA']]@data)==g),]
    # how many express above 0?
    # if more than 3, return TRUE, else return FALSE
    nexcells <- length(which(excells>0))
    
    if(nexcells>=3) {
      cat('\t \t \t \t \t \t >>> EXPRESSED \n')
      return(TRUE)
    } else {
      cat('\t \t \t \t \t \t >>> less than 3 cells expressing \n')
      return(FALSE)
    }
  })
  
})

# thno is a list, each slot is a dev stage, each element is TRUE/FALSE for presence/absence of the gene
# make it a dataframe
thno <- as.data.frame(do.call('rbind', thno))

# add dev stage as first column
thno <- thno %>%
  mutate(dev=names(seus), .before=1)

# change back zgc:172302 to abca7
colnames(thno)[which(colnames(thno)=='zgc:172302')] <- 'abca7'


### CHECKPOINT ####
# write.csv(thno, here('zfBrain_scRNAseq', 'thereorno.csv'), row.names=FALSE)
# import back
thno <- read.csv(here('zfBrain_scRNAseq', 'thereorno.csv'))

# any punctuation was swapped to .
# correct those gene names
colnames(thno)[which(colnames(thno)=='si.ch211.260p9.3')] <- 'si:ch211-260p9.3'
colnames(thno)[which(colnames(thno)=='zgc.174164')] <- 'zgc:174164'

# pivot longer for plot
thnol <- thno %>%
  pivot_longer(-dev,
               names_to='gene',
               values_to='there')

thnol$dev <- factor(thnol$dev, levels=rev(c('hpf12', 'hpf14', 'hpf16', 'hpf18', 'hpf20', 'hpf24', 'hpf36', 'dpf2', 'dpf3', 'dpf5', 'dpf8', 'dpf15')))

# make gene names italic
xlabs <- sapply(unique(thnol$gene), function(ge) {
  bquote(italic(.(ge)))
})

ggthno <- ggplot(thnol, aes(x=gene, y=dev)) +
  geom_tile(aes(fill=there), colour='white', linewidth=0.5) +
  scale_fill_manual(values=c('#d53e4f', '#abdda4')) +
  scale_x_discrete(position='top', labels=xlabs) +
  scale_y_discrete(labels=rev(c('12 hpf', '14 hpf', '16 hpf', '18 hpf', '20 hpf', '24 hpf', '36 hpf', '2 dpf', '3 dpf', '5 dpf', '8 dpf', '15 dpf'))) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x.top=element_text(angle=90, size=7, vjust=0.5, hjust=0, margin=margin(t=0, r=0, b=-2, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=-2, b=0, l=0)),
    legend.position='none',
  )
ggthno
ggsave(filename=here('zfBrain_scRNAseq', 'thereOrNo.pdf'), plot=ggthno, width=98, height=55, units='mm')

# label human genes?
# the only pairs with a different name are:
# ADAM10 : zgc:174164
# FCER1G : cd247l
# PLCG2 : si:ch211-260p9.3
# TSPOAP1 : CABZ01076737.1
# either write in legend, or can put in order of human genes on plot and add labels manually



# eLife reviews: how many genes total are expressed? ----------------------

# below takes a long time, see import of prepared file below

# minimum 3 cells (i.e. 4 or more) with expression
gencounts <- lapply(1:length(seus), function(st) {
  # loop through developmental stages
  cat('\n \t \t \t \t >>> Developmental stage:', names(seus)[st], '\n')
  
  tot <<- 0
  
  # for each developmental stage,
  # loop through genes present, and for each count how many cells express this gene > 0
  genpres <- row.names(seus[[st]][['RNA']]@data) # genes present
  
  sapply(1:length(genpres), function(g) {
    # count number of cells which express this gene
    cat('\t \t \t \t \t >>> gene', g, '/', length(genpres), '\n')
    gen <- genpres[g]
    # expression of every cell for gene g:
    excells <- seus[[st]][['RNA']]@data [which(rownames(seus[[st]][['RNA']]@data)==gen),]
    
    if(length(which(excells>0)) > 3) {
      tot <<- tot+1
    }
  })
  
  return(tot)
  
})

gencounts <- unlist(gencounts)
gencounts <- data.frame(names(seus), gencounts)
colnames(gencounts) <- c('devstage', 'ngenesexp')

### how many genes are there total?
# as maximum, will count any gene that is expressed at any timepoint
# loop through developmental stages again
allgenes <- lapply(1:length(seus), function(st) {
  # and for each, collect all the gene names
  genpres <- row.names(seus[[st]][['RNA']]@data) # genes present
  return(genpres)
})
# pool all the genes from all the developmental stages
allgenes <- unlist(allgenes)
# count unique genes
allgenes <- unique(allgenes)
length(allgenes) # 24,856

# add as column in gencounts
# then add ratio number of genes expressed / total number of genes that could be expressed
gencounts <- gencounts %>%
  mutate(ngenestot=length(allgenes)) %>%
  mutate(exp_pro=ngenesexp/ngenestot)

write.csv(gencounts, file=here('zfBrain_scRNAseq/nGenesExpressed.csv'), row.names=FALSE)


### re-import so we can skip time-consuming task above
gencounts <- read.csv(here('zfBrain_scRNAseq/nGenesExpressed.csv'))
# import AD genes present or not
thno <- read.csv(here('zfBrain_scRNAseq', 'thereorno.csv'))
# for each developmental stage, count how many AD genes are expressed out of total
adcounts <- sapply(1:nrow(thno), function(ro) {
  sum(thno[ro, 2:ncol(thno)])
})
# there are
ncol(thno)-1 # 42 AD genes
# add this info to gencounts
gencounts <- gencounts %>%
  mutate(nADexp=adcounts) %>%
  mutate(ntotAD=ncol(thno)-1) %>%
  mutate(expAD_pro=nADexp/ntotAD) %>%
  # add fold change
  mutate(foldAD=expAD_pro/exp_pro)

# make a barplot to illustrate
# to make stack barplots,
# need to add genes unexpressed proportions
gencounts <- gencounts %>%
  mutate(ngenesnotexp=ngenestot-ngenesexp, .after='ngenestot') %>%
  mutate(nADnotexp=ntotAD-nADexp, .after='ntotAD') %>%
  mutate(notexp_pro=1-exp_pro, .after='exp_pro') %>%
  mutate(notexpAD_pro=1-expAD_pro, .after='expAD_pro')

# simplify the data before pivot_longer
genco <- gencounts[,c('devstage', 'exp_pro', 'notexp_pro', 'expAD_pro', 'notexpAD_pro')]

colnames(genco) <- c('devstage', 'all_exp', 'all_notexp', 'AD_exp', 'AD_notexp')

gencol <- genco %>%
  pivot_longer(-devstage,
               names_to='set_exp',
               values_to='pro')
# now split set_exp to get all or AD sets
gencol <- gencol %>%
  mutate(set=strNthSplit(set_exp, '_', 1), .before='set_exp') %>%
  mutate(exp=strNthSplit(set_exp, '_', 2), .before='set_exp') %>%
  mutate(devstage_set=paste(devstage, set, sep='_'), .after='devstage')

gencol$set <- factor(gencol$set, levels=c('AD', 'all'))
gencol$exp <- factor(gencol$exp, levels=c('notexp', 'exp'))
gencol$set_exp <- factor(gencol$set_exp, levels=c('AD_notexp', 'AD_exp', 'all_notexp', 'all_exp'))
gencol$devstage <- factor(gencol$devstage, levels=c('hpf12', 'hpf14', 'hpf16', 'hpf18', 'hpf20', 'hpf24', 'hpf36',
                                                    'dpf2', 'dpf3', 'dpf5', 'dpf8', 'dpf15'))

expgenes <- ggplot(gencol, aes(x=set, y=pro, fill=set_exp)) +
  geom_col(width=0.9) +
  facet_grid(~devstage, scales='free_x', space='free') +
  scale_fill_manual(values=c('#d33d4e', '#abd3a3', '#e3828d', '#dfefdc')) +
  theme_minimal() +
  theme(
    strip.text.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(size=7, margin=margin(t=-5, r=0, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=-1, b=0, l=0)),
    legend.position='none') +
  coord_cartesian(ylim=c(0,1.0)) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels=c(0, 25, 50, 75, 100))
expgenes
ggsave(here('zfBrain_scRNAseq/expsets.pdf'), width=120, height=50, units='mm')

mean(gencounts$exp_pro)

# polished violin plots gene by gene --------------------------------------
# essentially reproduce violin plot from Seurat VlnPlot

# small function to prepare the data,
# then will plot manually gene by gene to choose colours etc.

# will do all 5 dpf as closest to HCR (6 dpf) and when behaviour experiment starts

prepViolin <- function(gene) {
  
  # get the data in a simple to use format
  gndf <- seu5dpf[['RNA']]@data [which(rownames(seu5dpf[['RNA']]@data)==gene),]
  gndf <- as.data.frame(gndf)
  # now 31659 rows, each is = one cell
  # 1 column
  # data = expression
  colnames(gndf) <- 'expr'
  
  # cell barcodes as rownames
  # put them as first column
  gndf <- gndf %>%
    mutate(cellid=row.names(.), .before=1)
  rownames(gndf) <- NULL
  
  # now add cluster assignment
  # from https://github.com/brlauuu/zf_brain/blob/master/data/stage_time_clustering_res_mapping.tsv
  # at 5 dpf, it is res.5.5
  # it is correct because there are 95 clusters (0--94)
  # other clusters have 53, 58, 63, 64, 68, 73, 81, 88, 101 clusters; so wrong
  # I think there is an error in cluster labels because cluster 94 (95th cluster) is not labelled
  # including in original file: https://github.com/brlauuu/zf_brain/blob/master/data/ZFBrainAtlasMaster.tsv
  # I think authors simply forgot to label it
  # will have to delete data from cluster 94 here
  clus5dpf <- seu5dpf$res.5.5
  
  # same numbers of cells in apoeb expression & cluster assignments
  if(length(clus5dpf) != nrow(gndf)) stop()
  
  # check the barcodes too
  if(!identical(sort(names(clus5dpf)), sort(gndf$cellid))) stop()
  
  # now, for each cell, find in which cluster it is
  gndf <- gndf %>%
    mutate(cluster=as.character(sapply(cellid, function(id) {
      clus5dpf[which(names(clus5dpf)==id)]
    })), .after='cellid')
  
  # add category to each cluster
  # which I wrote manually in zfclusters.csv
  # join all the columns to gndf, easier
  gndf <- left_join(gndf, cls[which(cls$stage=='5 dpf'),], by='cluster')
  
  # clusters currently in order 1, 10, 11, ...
  # sort them by number
  gndf$cluster <- factor(gndf$cluster, levels=mixedsort(unique(gndf$cluster)))
  # sort the dataframe this way
  gndf <- gndf[order(gndf$cluster),]
  
  # delete cells belonging to cluster #94
  # (see above for comments)
  gndf <- gndf[-which(gndf$cluster==94),]
  # reset the levels
  gndf$cluster <- factor(gndf$cluster)
  
  # set levels for the category
  # will worry about colours when plotting
  gndf$category <- factor(gndf$category)
  
  # return
  return(gndf)
  
}

#### apoeb ####

genexp <- prepViolin('apoeb')


# prepare vector of colours for categories
catcols <- c('#cb2a20', '#fcb505', '#588a9a', '#78ac63', '#982150', '#e47013', '#b2c4e0', '#697a87', '#db5072', '#417dcd')
names(catcols) <- levels(genexp$category)


# three clusters with high expression
# 40: epidermis (progenitors)
# 45: retina (muller glia)
# 58: microglia
# will colour them

# prepare colour vectors
# preallocate as all light grey
cols <- rep('#aeb3b4', length(unique(genexp$cluster)))
names(cols) <- sort(unique(genexp$cluster))
# and change the ones we want to colour
cols[which(names(cols)=='40')] <- '#fcb505'
cols[which(names(cols)=='45')] <- '#417dcd'
cols[which(names(cols)=='58')] <- '#cb2a20'

# now ready to plot
ggexp <- ggplot(genexp, aes(x=cluster, y=expr, colour=cluster)) +
  geom_quasirandom(width=0.4, size=0.6) +
  scale_colour_manual(values=cols) +
  theme_minimal() +
  theme(
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=5, angle=45, hjust=0, vjust=0, margin=margin(t=-2, r=0, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    legend.position='none'
  ) +
  ylab(expression(paste(italic('apoeb'), ' expression'))) +
  xlab('cluster')
ggexp
ggsave(here('zfBrain_scRNAseq', 'apoebExp.pdf'), width=162, height=50, units='mm', device=cairo_pdf)


#### appa ####

genexp <- prepViolin('appa')

# now ready to plot
ggexp <- ggplot(genexp, aes(x=cluster, y=expr, colour=category)) +
  geom_quasirandom(width=0.4, size=0.6) +
  scale_colour_manual(values=catcols) +
  theme_minimal() +
  theme(
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=5, angle=45, hjust=0, vjust=0, margin=margin(t=-2, r=0, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    legend.position='none'
  ) +
  ylab(expression(paste(italic('appa'), ' expression'))) +
  xlab('cluster')
ggexp
ggsave(here('zfBrain_scRNAseq', 'appaExp.pdf'), width=162, height=40, units='mm', device=cairo_pdf)


#### appb ####

genexp <- prepViolin('appb')

# now ready to plot
ggexp <- ggplot(genexp, aes(x=cluster, y=expr, colour=category)) +
  geom_quasirandom(width=0.4, size=0.6) +
  # scale_colour_manual(values=catcols) +
  theme_minimal() +
  theme(
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=5, angle=45, hjust=0, vjust=0, margin=margin(t=-2, r=0, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    legend.position='none'
  ) +
  ylab(expression(paste(italic('appb'), ' expression'))) +
  xlab('cluster')
ggexp
ggsave(here('zfBrain_scRNAseq', 'appbExp.pdf'), width=162, height=40, units='mm', device=cairo_pdf)


#### psen1 ####

genexp <- prepViolin('psen1')

# now ready to plot
ggexp <- ggplot(genexp, aes(x=cluster, y=expr, colour=category)) +
  geom_quasirandom(width=0.4, size=0.6) +
  scale_colour_manual(values=catcols) +
  theme_minimal() +
  theme(
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=5, angle=45, hjust=0, vjust=0, margin=margin(t=-2, r=0, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    legend.position='none'
  ) +
  ylab(expression(paste(italic('psen1'), ' expression'))) +
  xlab('cluster')
ggexp
ggsave(here('zfBrain_scRNAseq', 'psen1Exp.pdf'), width=162, height=40, units='mm', device=cairo_pdf)


#### psen2 ####

gen <- 'psen2'

genexp <- prepViolin(gen)

# now ready to plot
ggexp <- ggplot(genexp, aes(x=cluster, y=expr, colour=category)) +
  geom_quasirandom(width=0.4, size=0.6) +
  scale_colour_manual(values=catcols) +
  theme_minimal() +
  theme(
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=5, angle=45, hjust=0, vjust=0, margin=margin(t=-2, r=0, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    legend.position='none'
  ) +
  ylab(expression(paste(italic('psen2'), ' expression'))) +
  xlab('cluster')
ggexp

ggsave(here('zfBrain_scRNAseq', paste0(gen, 'Exp.pdf')), width=162, height=40, units='mm', device=cairo_pdf)


###### end of early onset
###### start of late onset


#### apoea ####

gen <- 'apoea'

genexp <- prepViolin(gen)

# now ready to plot
ggexp <- ggplot(genexp, aes(x=cluster, y=expr, colour=category)) +
  geom_quasirandom(width=0.4, size=0.6) +
  scale_colour_manual(values=catcols) +
  theme_minimal() +
  theme(
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=5, angle=45, hjust=0, vjust=0, margin=margin(t=-2, r=0, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    legend.position='none'
  ) +
  ylab(expression(paste(italic('apoea'), ' expression'))) +
  xlab('cluster')
ggexp

ggsave(here('zfBrain_scRNAseq', paste0(gen, 'Exp.pdf')), width=162, height=40, units='mm', device=cairo_pdf)


#### made apoeb above, will probably be in main figure


#### cd2ap ####

gen <- 'cd2ap'

genexp <- prepViolin(gen)

# now ready to plot
ggexp <- ggplot(genexp, aes(x=cluster, y=expr, colour=category)) +
  geom_quasirandom(width=0.4, size=0.6) +
  scale_colour_manual(values=catcols) +
  theme_minimal() +
  theme(
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=5, angle=45, hjust=0, vjust=0, margin=margin(t=-2, r=0, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    legend.position='none'
  ) +
  ylab(expression(paste(italic('cd2ap'), ' expression'))) +
  xlab('cluster')
ggexp

ggsave(here('zfBrain_scRNAseq', paste0(gen, 'Exp.pdf')), width=162, height=40, units='mm', device=cairo_pdf)


#### clu ####

gen <- 'clu'

genexp <- prepViolin(gen)

# now ready to plot
ggexp <- ggplot(genexp, aes(x=cluster, y=expr, colour=category)) +
  geom_quasirandom(width=0.4, size=0.6) +
  scale_colour_manual(values=catcols) +
  theme_minimal() +
  theme(
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=5, angle=45, hjust=0, vjust=0, margin=margin(t=-2, r=0, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    legend.position='none'
  ) +
  ylab(expression(paste(italic('clu'), ' expression'))) +
  xlab('cluster')
ggexp

ggsave(here('zfBrain_scRNAseq', paste0(gen, 'Exp.pdf')), width=162, height=40, units='mm', device=cairo_pdf)


#### sorl1 ####

gen <- 'sorl1'

genexp <- prepViolin(gen)

# now ready to plot
ggexp <- ggplot(genexp, aes(x=cluster, y=expr, colour=category)) +
  geom_quasirandom(width=0.4, size=0.6) +
  scale_colour_manual(values=catcols) +
  theme_minimal() +
  theme(
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=5, angle=45, hjust=0, vjust=0, margin=margin(t=-2, r=0, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    legend.position='none'
  ) +
  ylab(expression(paste(italic('sorl1'), ' expression'))) +
  xlab('cluster')
ggexp

ggsave(here('zfBrain_scRNAseq', paste0(gen, 'Exp.pdf')), width=162, height=40, units='mm', device=cairo_pdf)





# plot all genes/all stages -----------------------------------------------

thnol <- thno %>%
  pivot_longer(-gene,
               names_to='stage',
               values_to='there')


plotExpressionOneGeneAllStages <- function(gen) {

  Vlns <- vector(length=length(seus) , mode='list')
    
  # sum the booleans, TRUE = 1 / FALSE = 0, and gene name will turn into NA by as.logical (& removed by na.rm)
  # so will return number of developmental stage where expressed
  
  for(st in 1:length(seus)) { # loop through developmental stage
    # check that the gene is in the data at that stage, if not skip
    if (!as.logical(subset(thnol, gene==gen & stage==names(seus)[st], there))) { # if gene is *not* in the data
      cat('\t \t \t \t >>> gene', gen, 'is not in data at stage', names(seus)[st],'\n')
      
      Vlns[[st]] <- ggplot() + # add a blank canvas
        theme_void() +
        ggtitle(names(seus)[[st]]) +
        theme(
          plot.title=element_text(hjust=0.5, size=15)
        )
      
    } else {
      cat('\t \t \t \t >>> Plotting expression of gene', gen, 'at stage', names(seus)[st], '\n')
      Vlns[[st]] <- 
        VlnPlot(seus[[st]], features=gen) +
        ylab('expression') +
        ggtitle(names(seus)[st]) +
        theme(
          legend.position='none',
          axis.title.x=element_blank()
        )
    }
  }
  
  return(Vlns)
  
}


dir.create(here('zfBrain_scRNAseq', 'grids'))
# loop the function through the genes
for (g in 1:length(adgenes)) {
  exprGrid <- ggarrange(plotlist=plotExpressionOneGeneAllStages(adgenes[g]), nrow=length(seus), ncol=1)
  ggsave(filename=here('zfBrain_scRNAseq', 'grids', paste0(adgenes[g], '.pdf')), plot=exprGrid, width=400, height=1000, units='mm')
}



# try plotting with category colour ---------------------------------------

# just did 5 dpf (stage10) for now

# import file where I wrote that down
clscat <- read.csv(here('zfclusters_categories.tsv'),
                   sep='\t') # clusters

colnames(clscat) <- c('none', 'stage', 'cluster', 'match', 'celltype', 'category', 'markers')



# assign colours

apply(subset(clscat, stage=='5 dpf', category), 2, function(cat){
  
  if(cat=='nervous system') return('blue')
  if(cat=='else') return('grey')
  if(cat=='progenitors') return('green')
  if(cat=='eye') return('yellow')
  if(cat=='glia') return('orange')
  if(cat=='blood') return('red')
  
})


cats <- as.character(unlist(subset(clscat, stage=='5 dpf', category)))
catscol <- vector(mode='character', length=length(cats))

cat2col <- as.data.frame(matrix(nrow=length(unique(cats)), ncol=2))
colnames(cat2col) <- c('category', 'colour')

for (ca in 1:length(cats)) {
  if(ca=='nervous system')
  if(ca=='else') return('grey')
  if(ca=='progenitors') return('green')
  if(ca=='eye') return('yellow')
  if(ca=='glia') return('orange')
  if(ca=='blood') return('red')
}

VlnPlot(seus[[10]], features='apoeb') +
  scale_colour_manual(values=)
  
  ylab('expression') +
  ggtitle(names(seus)[10]) +
  scale_colour_manual() +
  theme(
    legend.position='none',
    axis.title.x=element_blank())


allgns <- rownames(seus[[10]])
pbmc <- ScaleData(seus[[10]], features=allgns, vars=scale)

dim(seus[[10]][['RNA']]@scale.data)



# plot the tSNE -----------------------------------------------------------

FeaturePlot(seus[[1]], features='apoeb', label=TRUE)