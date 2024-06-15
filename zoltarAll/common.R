#####################################################
# ~ ZFAD: explore possible overlaps between ZOLTAR predictions ~
#
# for eLife reviews
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################


# packages ----------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(tibble)

library(ggplot2)


# import ------------------------------------------------------------------

genes <- c('appab', 'psen1', 'psen2', 'apoeab', 'cd2ap', 'clu', 'sorl1')


## import drugs ranked in a list
drul <- lapply(genes, function(gen) {
  
  filnm <- paste0(gen, '_drugsRanked.csv')
  return( read.csv(here('zoltarAll', filnm)) )
  
})
names(drul) <- genes


## import indications in a list
indl <- lapply(genes, function(gen) {
  
  filnm <- paste0(gen, '_indications.csv')
  return( read.csv(here('zoltarAll', filnm)) )
  
})
names(indl) <- genes


## import targets in a list
tarl <- lapply(genes, function(gen) {
  
  filnm <- paste0(gen, '_TTDtargets.csv')
  return( read.csv(here('zoltarAll', filnm)) )
  
})
names(tarl) <- genes


## import KEGG pathways in a list
kegl <- lapply(genes, function(gen) {
  
  filnm <- paste0(gen, '_KEGGpathways.csv')
  return( read.csv(here('zoltarAll', filnm)) )
  
})
names(kegl) <- genes



# indications -------------------------------------------------------------

# for each indication,
# count how many genes have it as significant & record them
uano <- indl[[1]]$annotation
# uano for unique annotation

# tal for tally
tal <- sapply(1:length(uano), function(ui) {
  
  # annotation is:
  ano <- uano[ui]
  
  unlist(lapply(1:length(indl), function(g) {
    
    # get that genes's data
    gd <- indl[[g]]
    
    # find that annotation (the row index)
    ani <- which(gd$annotation==ano)
    
    # if pval is below 0.05, return the gene name
    if(!is.na(gd[ani, 'pval']) & gd[ani, 'pval'] < 0.05) {
      return(names(indl)[g])
    }
    
  }))
  
})
names(tal) <- uano

# we now have a list, each slot is an indication
# each element is a vector with the genes that have this annotation

# now, simply count annotation by annotation number of genes
# cos for counts

cos <- unlist(lapply(tal, length))
cos <- data.frame(annotation=names(cos), count=cos)
row.names(cos) <- NULL

# order from largest count to smallest
cos <- cos[rev(order(cos$count)),]

# can do a small plot
cos$annotation <- factor(cos$annotation, levels=cos$annotation)

# skip all the 0
cosNo0 <- cos %>%
  filter(count>0)

ggplot(cosNo0, aes(x=annotation, y=count)) +
  geom_bar(stat='identity') +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title.x=element_blank(),
    panel.grid.minor.y=element_blank()
  ) +
  coord_cartesian(ylim=c(0, 7)) +
  scale_y_continuous(breaks=0:7)

# ideal would be to have colours with the gene names
# for this, we need data in long format that says e.g.:
# insomnia / cd2ap / 0
# insomnia / clu / 1
# etc.
# tal is already pretty close

# cannot find a straightforward solution...

# go through tal one annotation at a time
# and prepare row of final dataframe
# which should be: one column per gene, one row per annotation
# 1 if significant, 0 if not
copg <- lapply(tal, function(ta) {
  sapply(genes, function(gen) {
    if(gen %in% ta) return(1)
    if(!gen %in% ta) return(0)
  })
})
# copg for count per gene
# now rbind into a dataframe
copg <- as.data.frame(do.call(rbind, copg))
# do not use row.names
copg <- copg %>%
  mutate(annotation=row.names(copg), .before=1)
row.names(copg) <- NULL

# the levels of annotation need to be most frequent to least frequent
# we can simply copy the order from cos
copg$annotation <- factor(copg$annotation, levels=cos$annotation)

# remove the annotation that are 0
# can just copy from cos as well
copgNo0 <- copg[copg$annotation %in% cosNo0$annotation,]

## pivot long
copgl <- copgNo0 %>%
  pivot_longer(-annotation,
               names_to='gene',
               values_to='count')

## plot
ggCo <- ggplot(copgl, aes(x=annotation, y=count, fill=gene)) +
  geom_bar(stat='identity') +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title.x=element_blank(),
    panel.grid.minor.y=element_blank()
  ) +
  coord_cartesian(ylim=c(0, 7)) +
  scale_y_continuous(breaks=0:7) +
  ylab('number of genes for which significant')
ggsave(here('zoltarAll/plots/indicationCount.pdf'), ggCo, width=200, height=200, units='mm')


### maybe better is heatmap

# reverse the order of annotations,
# I want most frequent on top
copgl$annotation <- factor(copgl$annotation, levels=rev(levels(copgNo0$annotation)))
# order of genes:
copgl$gene <- factor(copgl$gene, levels=c('appab', 'psen1', 'psen2',
                                          'apoeab', 'cd2ap', 'clu', 'sorl1'))

copgl$count <- factor(copgl$count)

ggHea <- ggplot(copgl, aes(x=gene, y=annotation, fill=count)) +
  geom_tile(colour='white', linewidth=0.5) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(angle=90, size=7, hjust=0),
    axis.text.y=element_text(size=7),
    legend.title=element_blank(),
    legend.text=element_text(size=5),
    plot.margin=unit(c(0, 0, 0, 0), 'pt'),
    legend.position='none'
  ) +
  scale_x_discrete(position='top') +
  coord_equal() +
  scale_fill_manual(values=c('#aeb3b4', '#cb2a20'))
ggHea
ggsave(here('zoltarAll/plots/indicationHeat.pdf'), ggHea, width=100, height=200, units='mm')


# targets -----------------------------------------------------------------
# just repeating the code above

# for each target,
# count how many genes have it as significant & record them
uano <- tarl[[1]]$TARGNAME
# uano for unique annotation

# tal for tally
tal <- sapply(1:length(uano), function(ui) {
  
  # annotation is:
  ano <- uano[ui]
  
  unlist(lapply(1:length(tarl), function(g) {
    
    # get that genes's data
    gd <- tarl[[g]]
    
    # find that annotation (the row index)
    ani <- which(gd$TARGNAME==ano)
    
    # if pval is below 0.05, return the gene name
    if(!is.na(gd[ani, 'pval']) & gd[ani, 'pval'] < 0.05) {
      return(names(tarl)[g])
    }
    
  }))
  
})
names(tal) <- uano

# we now have a list, each slot is an annotation
# each element is a vector with the genes that have this annotation

# now, simply count annotation by annotation number of genes
# cos for counts

cos <- unlist(lapply(tal, length))
cos <- data.frame(annotation=names(cos), count=cos)
row.names(cos) <- NULL

# order from largest count to smallest
cos <- cos[rev(order(cos$count)),]

# can do a small plot
cos$annotation <- factor(cos$annotation, levels=cos$annotation)

# skip all the 0
cosNo0 <- cos %>%
  filter(count>0)

ggplot(cosNo0, aes(x=annotation, y=count)) +
  geom_bar(stat='identity') +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title.x=element_blank(),
    panel.grid.minor.y=element_blank()
  ) +
  coord_cartesian(ylim=c(0, 7)) +
  scale_y_continuous(breaks=0:7)

# ideal would be to have colours with the gene names
# for this, we need data in long format that says e.g.:
# insomnia / cd2ap / 0
# insomnia / clu / 1
# etc.
# tal is already pretty close

# cannot find a straightforward solution...

# go through tal one annotation at a time
# and prepare row of final dataframe
# which should be: one column per gene, one row per annotation
# 1 if significant, 0 if not
copg <- lapply(tal, function(ta) {
  sapply(genes, function(gen) {
    if(gen %in% ta) return(1)
    if(!gen %in% ta) return(0)
  })
})
# copg for count per gene
# now rbind into a dataframe
copg <- as.data.frame(do.call(rbind, copg))
# do not use row.names
copg <- copg %>%
  mutate(annotation=row.names(copg), .before=1)
row.names(copg) <- NULL

# the levels of annotation need to be most frequent to least frequent
# we can simply copy the order from cos
copg$annotation <- factor(copg$annotation, levels=cos$annotation)

# remove the annotation that are 0
# can just copy from cos as well
copgNo0 <- copg[copg$annotation %in% cosNo0$annotation,]

## pivot long
copgl <- copgNo0 %>%
  pivot_longer(-annotation,
               names_to='gene',
               values_to='count')

## plot
ggCo <- ggplot(copgl, aes(x=annotation, y=count, fill=gene)) +
  geom_bar(stat='identity') +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title.x=element_blank(),
    panel.grid.minor.y=element_blank()
  ) +
  coord_cartesian(ylim=c(0, 7)) +
  scale_y_continuous(breaks=0:7) +
  ylab('number of genes for which significant')
ggCo
ggsave(here('zoltarAll/plots/targetCount.pdf'), ggCo, width=200, height=200, units='mm')


### maybe better is heatmap

# reverse the order of annotations,
# I want most frequent on top
copgl$annotation <- factor(copgl$annotation, levels=rev(levels(copgNo0$annotation)))
# order of genes:
copgl$gene <- factor(copgl$gene, levels=c('appab', 'psen1', 'psen2',
                                          'apoeab', 'cd2ap', 'clu', 'sorl1'))

copgl$count <- factor(copgl$count)

ggHea <- ggplot(copgl, aes(x=gene, y=annotation, fill=count)) +
  geom_tile(colour='white', linewidth=0.5) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(angle=90, size=7, hjust=0),
    axis.text.y=element_text(size=7),
    legend.title=element_blank(),
    legend.text=element_text(size=5),
    plot.margin=unit(c(0, 0, 0, 0), 'pt'),
    legend.position='none'
  ) +
  scale_x_discrete(position='top') +
  coord_equal() +
  scale_fill_manual(values=c('#aeb3b4', '#cb2a20'))
ggHea
ggsave(here('zoltarAll/plots/targetHeat.pdf'), ggHea, width=100, height=200, units='mm')



# KEGG pathways -----------------------------------------------------------

# just repeating the code above

# for each KEGG,
# count how many genes have it as significant & record them
# issue here is KEGG name has a few duplicates
kegl[[1]][which(duplicated(kegl[[1]]$keggname)),]

# I think duplicates are KEGG pathways for other species than humans
# e.g.:
# hsa is Homo sapiens
# pfa is Plasmodium falciparum
# sce is Saccharomyces cerevisiae
# etc
# from https://www.genome.jp/kegg-bin/show_organism?orgs=sce

# will delete manually here, but should do it in ZOLTAR app

uano <- kegl[[1]]$annotation
# uano for unique annotation

# tal for tally
tal <- sapply(1:length(uano), function(ui) {
  
  # annotation is:
  ano <- uano[ui]
  
  unlist(lapply(1:length(kegl), function(g) {
    
    # get that genes's data
    gd <- kegl[[g]]
    
    # find that annotation (the row index)
    ani <- which(gd$annotation==ano)
    
    # if pval is below 0.05, return the gene name
    if(!is.na(gd[ani, 'pval']) & gd[ani, 'pval'] < 0.05) {
      return(names(kegl)[g])
    }
    
  }))
  
})
names(tal) <- uano


# we now have a list, each slot is an annotation
# each element is a vector with the genes that have this annotation

# now, simply count annotation by annotation number of genes
# cos for counts

cos <- unlist(lapply(tal, length))
cos <- data.frame(annotation=names(cos), count=cos)
row.names(cos) <- NULL

# order from largest count to smallest
cos <- cos[rev(order(cos$count)),]

# add back the name
# we could not use it earlier because not unique
# e.g. there are two "Biosynthesis of antibiotics", for some reason
kegnms <- sapply(cos$annotation, function(ano) {
  return( kegl[[1]][ which(kegl[[1]]$annotation==ano) , 'keggname'] )
})
cos <- cos %>%
  mutate(keggname=kegnms, .after='annotation')

# DELETE OTHER THAN HUMAN
cos <- cos[startsWith(cos$annotation, 'hsa'),]

# now should be unique annotations
sum(duplicated(cos$annotation))

# can do a small plot
cos$keggname <- factor(cos$keggname, levels=cos$keggname)


# skip all the 0
cosNo0 <- cos %>%
  filter(count>0)

ggplot(cosNo0, aes(x=keggname, y=count)) +
  geom_bar(stat='identity') +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title.x=element_blank(),
    panel.grid.minor.y=element_blank()
  ) +
  coord_cartesian(ylim=c(0, 7)) +
  scale_y_continuous(breaks=0:7)

# ideal would be to have colours with the gene names
# for this, we need data in long format that says e.g.:
# insomnia / cd2ap / 0
# insomnia / clu / 1
# etc.
# tal is already pretty close

# go through tal one annotation at a time
# and prepare row of final dataframe
# which should be: one column per gene, one row per annotation
# 1 if significant, 0 if not
copg <- lapply(tal, function(ta) {
  sapply(genes, function(gen) {
    if(gen %in% ta) return(1)
    if(!gen %in% ta) return(0)
  })
})
# copg for count per gene
# now rbind into a dataframe
copg <- as.data.frame(do.call(rbind, copg))
# do not use row.names
copg <- copg %>%
  mutate(annotation=row.names(copg), .before=1)
row.names(copg) <- NULL

# add back the name
# we could not use it earlier because not unique
# e.g. there are two "Biosynthesis of antibiotics", for some reason
kegnms <- sapply(copg$annotation, function(ano) {
  return( kegl[[1]][ which(kegl[[1]]$annotation==ano) , 'keggname'] )
})
copg <- copg %>%
  mutate(keggname=kegnms, .after='annotation')

# DELETE OTHER THAN HUMAN
copg <- copg[startsWith(copg$annotation, 'hsa'),]

# now should be unique annotations
sum(duplicated(copg$annotation))

# the levels of annotation need to be most frequent to least frequent
# we can simply copy the order from cos
copg$annotation <- factor(copg$annotation, levels=cos$annotation)
copg$keggname <- factor(copg$keggname, levels=cos$keggname)

# remove the annotation that are 0
# can just copy from cos as well
copgNo0 <- copg[copg$annotation %in% cosNo0$annotation,]

## pivot long
copgl <- copgNo0 %>%
  pivot_longer(-c(annotation, keggname),
               names_to='gene',
               values_to='count')

## plot
ggCo <- ggplot(copgl, aes(x=keggname, y=count, fill=gene)) +
  geom_bar(stat='identity') +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title.x=element_blank(),
    panel.grid.minor.y=element_blank()
  ) +
  coord_cartesian(ylim=c(0, 7)) +
  scale_y_continuous(breaks=0:7) +
  ylab('number of genes for which significant')
ggCo
ggsave(here('zoltarAll/plots/KEGGCount.pdf'), ggCo, width=200, height=200, units='mm')


### maybe better is heatmap

# reverse the order of annotations,
# I want most frequent on top
copgl$keggname <- factor(copgl$keggname, levels=rev(levels(copgNo0$keggname)))
# order of genes:
copgl$gene <- factor(copgl$gene, levels=c('appab', 'psen1', 'psen2',
                                          'apoeab', 'cd2ap', 'clu', 'sorl1'))

copgl$count <- factor(copgl$count)

ggHea <- ggplot(copgl, aes(x=gene, y=keggname, fill=count)) +
  geom_tile(colour='white', linewidth=0.5) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(angle=90, size=7, hjust=0),
    axis.text.y=element_text(size=7),
    legend.title=element_blank(),
    legend.text=element_text(size=5),
    plot.margin=unit(c(0, 0, 0, 0), 'pt'),
    legend.position='none'
  ) +
  scale_x_discrete(position='top') +
  coord_equal() +
  scale_fill_manual(values=c('#aeb3b4', '#cb2a20'))
ggHea
ggsave(here('zoltarAll/plots/KEGGHeat.pdf'), ggHea, width=100, height=200, units='mm')
