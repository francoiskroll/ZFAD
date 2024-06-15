#####################################################
# ~ ZFAD: create orthologues donut plot ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################



# packages ----------------------------------------------------------------

library(here)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(tidyr)



# import ------------------------------------------------------------------

# import database written by buildZFADdb.R
ort <- read.xlsx(here('geneSelect', 'ZFAD.xlsx'), sheet='zebrafish_orthologues')



# count how many orthologues for each human gene --------------------------

# ! if proceed to count how many orthologues per human gene below, will count NA as one, but should be zero
# first exclude the NA
ortnoNA <- ort[!is.na(ort$zebrafishgene),]

# then count how many orthologues per human gene
ortn <- ortnoNA %>%
  group_by(humangene) %>%
  summarise_at(vars(zebrafishgene),
               list(
                 north= ~ length(.)
               ))

# and now add the count 0
# human genes without a zebrafish orthologue:
h_noort <- ort[which(is.na(ort$zebrafishgene)), 'humangene']
ortn <- rbind(ortn, data.frame(humangene=h_noort, north=0))

# check we did not create any duplicates
sum(duplicated(ortn$humangene))

# check all the human genes are there
identical(sort(unique(ort$humangene)), sort(ortn$humangene))

# OK, good to go

# add column,
# if 1: one2one
# if 2: one2two
# if 3+: one2many
ortn <- ortn %>%
  mutate(type=sapply(north, function(no){
    if(no==0) return('one2none')
    if(no==1) return('one2one')
    if(no==2) return('one2two')
    if(no>=3) return('one2many')
  }))

# now count how many humangenes in each category
# i.e. how many one2one, etc.
# co for counts
ortco <- ortn %>%
  group_by(type) %>%
  tally()

# goal is to prepare data as if making a stacked barplot
# then we make the grid a circle, essentially
# I am following https://r-graph-gallery.com/128-ring-or-donut-plot.html

# compute proportions
ortco$pro = ortco$n / sum(ortco$n)

# compute the cumulative proportions
# (top of each section in the stacked barplot)
ortco$ymax = cumsum(ortco$pro)

# compute the bottom of each section
ortco$ymin = c(0, head(ortco$ymax, n=-1))

# set the levels
ortco$type <- factor(ortco$type, levels=c('one2none', 'one2one', 'one2two', 'one2many'))

# draw the plot
orthdonut <- ggplot(ortco, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
  geom_rect(colour='black') +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void() +
  scale_fill_manual(values=c('#d7dbe1', '#5f717e', '#859070', '#ec9278')) +
  theme(legend.position='none')
orthdonut

ggsave(here('geneSelect', 'orthoDonut.pdf'), orthdonut, width=47, height=47, units='mm', device=cairo_pdf)


### make slices
# add a column that is all 1, as each gene should contribute same slice to the pie
ortn$slice <- 1
# compute proportions
ortn$pro = ortn$slice / sum(ortn$north)

# set the levels
ortn$type <- factor(ortn$type, levels=c('one2none', 'one2one', 'one2two', 'one2many'))

# sort data according to levels
ortn <- ortn[order(ortn$type),]

# compute the cumulative proportions
# (top of each section in the stacked barplot)
ortn$ymax = cumsum(ortn$pro)

# compute the bottom of each section
ortn$ymin = c(0, head(ortn$ymax, n=-1))

# set genes as factor in this order, otherwise it plots alphabetically
ortn$humangene <- factor(ortn$humangene, levels=ortn$humangene)

# draw the plot
orthdonut <- ggplot(ortn, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
  # geom_rect(colour='black') +
  geom_rect(colour='#eaedf0', linewidth=0.1) +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position='none') +
  scale_fill_manual(values=c('#c3cad3', '#5f717e', '#859070', '#ec9278'))
orthdonut

ggsave(here('geneSelect', 'orthoDonutv2.pdf'), orthdonut, width=47, height=47, units='mm', device=cairo_pdf)


### how many human genes?
sum(ortco$n)
nrow(ortn)
