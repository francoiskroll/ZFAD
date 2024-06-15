# transcriptMapper.R
# plots a transcript schematics from any Ensembl transcript ID
# zebrafish only, but would just need to change useDataset(...) for another species

# v1

# v2
# made it possible to not mark any position by giving positions_tomark=NA
# added various checks/error messages, e.g. to check it actually found the transcript ID in Ensembl or to check it found positions to mark
# allowed no left UTR and/or no right UTR. I had never seen the case before so did not think about coding it

# v3
# moved to utilities/ folder of ZFAD

# v4
# PAM could be exactly first or last nucleotide of exon, did not think about that before
# v3 and before was always saying it was finding all positions because it was not clearing 'mrk' from finding previous position
# consequently, if e.g. it could not find third (of 3) positions, it would return position1, position2, position2

# v5
# some minor corrections to allow transcript with no UTR whatsoever
# I had never seen this case before
# if need example: ENSDART00000166786


# packages & small functions ----------------------------------------------

# for installations, you might need (in that order)

# install.packages('BiocManager')
# BiocManager::install("biomaRt")

# install.packages('devtools')
# devtools::install_github("hrbrmstr/statebins")

library(biomaRt)
library(ggplot2)
library(dplyr)
library(tidyr)



# function transcriptMapper(...) ------------------------------------------

# transcript_id = Ensembl transcript ID, e.g. 'ENSDART00000172219'

# positions_tomark = any one or multiple genomic position to mark on the map, e.g. c(51218545, 51216830, 51210589)
# give NA to not mark any position, as such: positions_tomark = NA

# utr_height = height of 5' and 3'-UTR (arbitrary units)
# exon_height = height of protein-coding exons (arbitrary units)
# intron_height = height of introns (arbitrary units)
# you want exon_height > utr_height > intron_height; e.g. 8 > 5 > 2

# scale_exons = by how much should we scale down exon lengths, e.g. * 0.1, i.e. exon lengths will be a tenth of the real length
# scale_introns = by how much should we scale down intron lengths, e.g. * 0.01, i.e. exon lengths will be a hundredth of the real length

# exon_colour = fill colour of exons
# intron_colour = fill colour of introns

# roundCorner = TRUE or FALSE, whether exons should have round corners
# corner_radius = radius of that corner (only applies if roundCorner=TRUE, may need to do by trial and error when exporting in pdf)

# arrowsize = size of the arrows that mark the positions_tomark

# exportFolder = path to the export folder, by default filename will be transcript_id.pdf, e.g. ENSDART00000172219.pdf

# example

# transcriptMapper(transcript_id='ENSDART00000149864',
#                  
#                  positions_tomark=c(51218545, 51216830, 51210589),
#                  
#                  utr_height=5,
#                  exon_height=8,
#                  intron_height=2,
#                  
#                  scale_exons=0.1,
#                  scale_introns=0.01,
#                  
#                  exon_colour='#697a87',
#                  intron_colour='#595E60',
#                  
#                  roundCorners=TRUE,
#                  corner_radius=1.5,
#                  
#                  arrowsize=20,
#                  
#                  exportFolder=here())

transcriptMapper <- function(transcript_id,
                             
                             positions_tomark,
                             
                             utr_height=5,
                             exon_height=8,
                             intron_height=2,
                             
                             scale_exons=0.1,
                             scale_introns=0.01,
                             
                             exon_colour='#697a87',
                             intron_colour='#595E60',
                             
                             roundCorners=TRUE,
                             corner_radius=3,
                             
                             arrowsize=10,
                             arrowpos='above',
                             
                             exportFolder,
                             width=90,
                             height=15) {
  

  # warn user about Ensembl -------------------------------------------------

  cat('\n \t \t \t \t Disclaimer: accessing Ensembl is frequently glitchy. \
      \t \t \t \t If it stops running early and the error message looks like it might be from Ensembl, try running again a few times. \n \n')
  
  
  # get positions of UTR/exon/intron boundaries -----------------------------
  
  
  cat('\t \t \t \t >>> Searching Ensembl for intron/exon positions... \n')
  
  # set up the database
  mart <- useMart('ensembl')
  mart <- useDataset('drerio_gene_ensembl', mart)
  
  # 
  attr <- c('ensembl_transcript_id',
            'rank',
            '5_utr_start',
            '5_utr_end',
            'exon_chrom_start',
            'exon_chrom_end',
            '3_utr_start',
            '3_utr_end')
  
  tpo <- getBM(attributes=attr, filters='ensembl_transcript_id', values=transcript_id, mart=mart)
  
  # check we found the transcript
  if (nrow(tpo)==0) stop ('\n \n \t \t \t \t >>> Error: could not find this transcript ID in Ensembl. Check the following: \
    \t \t \t \t -- Are you using a transcript ID? It should start with ENSDART... \
    \t \t \t \t -- Do not include the version of the transcript ID. e.g. for ENSDART00000123174.4, only give ENSDART00000123174 \n')
  
  
  

  # order the exons ---------------------------------------------------------

  # for reasons I do not understand, Ensembl sometimes give the exons in a shuffled order
  # column `rank` scan be used to put them in order
  tpo <- tpo[order(tpo$rank),]
  
  
  # flip ordering if negative strand ----------------------------------------

  # is transcript + or -?
  if(tpo$exon_chrom_start[1] > tpo$exon_chrom_end[nrow(tpo)]) {positive <- FALSE} # if positions decrease, then negative strand
  if(tpo$exon_chrom_start[1] < tpo$exon_chrom_end[nrow(tpo)]) {positive <- TRUE} # if positions increase, then positive strand
  
  # if negative strand, flip order of table
  if(!positive) {
    
    tpo <- tpo[nrow(tpo):1 ,]
    
  }
  
  # in either case, make new rank column
  # keep original rank table
  colnames(tpo) [which(colnames(tpo)=='rank')] <- 'original_rank'
  
  # make new one
  tpo$rank <- 1:nrow(tpo)
  
  # from now on, will call left UTR and right UTR
  # if positive strand, left UTR = 5'-UTR; right UTR = 3'-UTR
  # if negative strand, left UTR = 3'UTR; right UTR = 5'-UTR
  # (as we are always drawing in 5'-3' genome direction)
  
  
  # format positions table --------------------------------------------------
  
  # preallocate new table
  # number of rows should be number of left UTRs + number of exons + number of introns + number of right UTRs
  
  # count UTRs
  # ! not always a left & right UTR, e.g. ENSDART00000188300 does not have a left UTR
  # ! can be splitted between two or more exons, which is counted as multiple UTRs here
  
  # how many left UTRs?
  if (positive) { nlutr <- sum(!is.na(tpo$`5_utr_start`)) }
  if (!positive) { nlutr <- sum(!is.na(tpo$`3_utr_start`)) }
  
  # how many right UTRs?
  if (positive) { nrutr <- sum(!is.na(tpo$`3_utr_start`)) }
  if (!positive) { nrutr <- sum(!is.na(tpo$`5_utr_start`)) }
  
  # how many exons?
  nex <- max(tpo$rank)
  
  # how many protein-coding exons?
  
  # there are four cases here
  
  # 1-- (most often) at least 1 left UTR & at least 1 right UTR
  if (nlutr>0 & nrutr>0) {
    nexpc <- nex - (nlutr-1) - (nrutr-1)
    # goal here is to only count protein-coding exons
    # to do so, we subtract non-protein-coding exons (i.e. full UTR exons) from the exon count
    # But, if an exon is part UTR/part protein-coding, we want to count as 1
    # hence nUTR minus 1, so we do not subtract the exon that is part UTR/part protein-coding ***
  
  # 2-- 0 left UTR & at least 1 right UTR, e.g. ENSDART00000188300  
  } else if (nlutr==0 & nrutr>0) {
    nexpc <- nex - (nrutr-1)
    # avoids to do nex - (0-1), which ends up incorrectly adding 1 protein-coding exon
    
  # 3-- at least 1 left UTR & 0 right UTR, I have never seen but coded it just in case
  } else if (nlutr>0 & nrutr==0) {
    nexpc <- nex - (nlutr-1)
    # avoids to do nex - (0-1), which ends up incorrectly adding 1 protein-coding exon
  
  # 4-- 0 left UTR & right UTR, I have never seen but will code it just in case
  } else if (nlutr==0 & nrutr==0) {
    nexpc <- nex
    # if no UTR, then all exons are protein-coding
  }
  
  # *** this does not allow a case where full UTR exon is followed by full coding exon
  # i.e. a 'clean' break between UTR and coding
  # I have never seen an example so will not code for now but might change in the future
  
  # first make a table with just UTRs and exons
  tpo2 <- as.data.frame(matrix(nrow= nlutr + nexpc + nrutr, ncol=7))
  
  colnames(tpo2) <- c('type', 'feature', 'feature_num', 'start', 'stop', 'ystart', 'ystop')
  
  # fill in feature
  feat <- c( rep('leftutr', nlutr) ,
             rep('exon', nexpc) ,
             rep('rightutr', nrutr))
  tpo2$feature <- feat
  
  # where do we need to insert introns?
  # after every row
  inin <- 1:nrow(tpo2) # inin = insert introns (after these row indices)
  
  # except if left UTR then exon or exon then right UTR
  exceptrow <- c()
  
  # check every exon for this case
  for (i in 1:length(feat)) {
    
    # ! careful to not go after last feature otherwise it will throw an error
    # this can happen if there is no UTR whatsoever
    if(i+1 > length(feat)) next()
    
    # do we have left UTR followed by exon?
    if(feat[i]=='leftutr' & feat[i+1]=='exon') {exceptrow <- c(exceptrow, i)}
    # do we have exon followed by right UTR?
    if(feat[i]=='exon' & feat[i+1]=='rightutr') {exceptrow <- c(exceptrow, i)}
  }
  # remove these indices, if any
  # i.e. we do not want to add an intron after those features
  if(length(exceptrow)>0) {
    inin <- inin[-which(inin %in% exceptrow)]
  }

  # ! we should also not add an intron after the last exon
  # i.e. transcript finishes with an exon or UTR, not an intron
  inin <- inin[-length(inin)]
  
  # insert intron rows
  for (i in 1:length(inin)) {
    tpo2 <- tpo2 %>%
      add_row(feature='intron', .after=inin[i])
    
    inin <- inin + 1 # need to update the indices each time as we are shifting everything
  }
  
  # fill in type
  tpo2$type <- sapply(tpo2$feature, function(f) {
    if (f=='leftutr') return ('exon')
    if (f=='rightutr') return ('exon')
    if (f=='exon') return ('exon')
    if (f=='intron') return ('intron')
  })
  
  # now replace any exon in feature so it says proteincoding exon, abbreviated as pcexon
  tpo2[which(tpo2$feature=='exon'), 'feature'] <- 'pcexon'
  
  # ystart (bottom line) and ystop depend on the feature
  # will make line centered on 10
  
  # for left UTR and right UTR, will go from 10 - (utr_height/2) up to 10 + (utr_height/2)
  tpo2[which(tpo2$feature=='leftutr' | tpo2$feature=='rightutr') , 'ystart'] <- 10 - (utr_height/2)
  tpo2[which(tpo2$feature=='leftutr' | tpo2$feature=='rightutr') , 'ystop'] <- 10 + (utr_height/2)
  
  # for exon (but not UTR), will go from 10 - (exon_height/2) up to 10 + (exon_height/2)
  tpo2[which(tpo2$feature=='pcexon') , 'ystart'] <- 10 - (exon_height/2)
  tpo2[which(tpo2$feature=='pcexon') , 'ystop'] <- 10 + (exon_height/2)
  
  # for intron, will go from 10 - (intron_height/2) up to 10 + (intron_height/2)
  tpo2[which(tpo2$type=='intron') , 'ystart'] <- 10 - (intron_height/2)
  tpo2[which(tpo2$type=='intron') , 'ystop'] <- 10 + (intron_height/2)
  
  
  # now we will fill start & stop positions

  # start with 5'-UTR & 3'-UTR
  
  # take 5'-UTR start, 5'-UTR stop, 3'-UTR start, 3'-UTR stop
  
  if (positive) {
    
    lutr_start <- tpo$`5_utr_start`[!is.na(tpo$`5_utr_start`)]
    lutr_stop <- tpo$`5_utr_end`[!is.na(tpo$`5_utr_end`)]
    rutr_start <- tpo$`3_utr_start`[!is.na(tpo$`3_utr_start`)]
    rutr_stop <- tpo$`3_utr_end`[!is.na(tpo$`3_utr_end`)]
    
  } else if (!positive) {
    
    lutr_start <- tpo$`3_utr_start`[!is.na(tpo$`3_utr_start`)]
    lutr_stop <- tpo$`3_utr_end`[!is.na(tpo$`3_utr_end`)]
    rutr_start <- tpo$`5_utr_start`[!is.na(tpo$`5_utr_start`)]
    rutr_stop <- tpo$`5_utr_end`[!is.na(tpo$`5_utr_end`)]
    
  }
  
  # we place them in tpo2
  tpo2[which(tpo2$feature=='leftutr'), 'start'] <- lutr_start
  tpo2[which(tpo2$feature=='leftutr'), 'stop'] <- lutr_stop
  tpo2[which(tpo2$feature=='rightutr'), 'start'] <- rutr_start
  tpo2[which(tpo2$feature=='rightutr'), 'stop'] <- rutr_stop
  
  
  # now for exons
  
  # will be easier to first correct original dataframe, so that
  # the start position of the first protein-coding exon is the end of the last 5'-UTR
  # the end position of exon n should be the start of the first 3'-UTR
  
  # before that,
  # are there more than one left UTR? then first n - 1 rows are entirely UTR exons, should delete them from original dataframe
  # are there more than one right UTR? then last n - 1 rows are entirely UTR exons, should delete them from original dataframe
  if (nlutr > 1) {tpo <- tpo[- head(1:nrow(tpo), nlutr-1) , ]}
  if (nrutr > 1) {tpo <- tpo[- tail(1:nrow(tpo), nrutr-1) , ]}
  
  # now correct as above
  
  # ! there may be no left UTR, in which case we do not change the start position of the first protein-coding exon
  if (nlutr > 0) { tpo[which(tpo$rank==1), 'exon_chrom_start'] <- lutr_stop[nlutr] }
  
  # ! there may be no right UTR, in which case we do not change the end position of the last protein-coding exon
  if (nrutr > 0) { tpo[which.max(tpo$rank), 'exon_chrom_end'] <- rutr_start[1] }
  
  # now we place the exon positions in tpo2
  tpo2[which(tpo2$feature=='pcexon') , 'start'] <- tpo$exon_chrom_start
  tpo2[which(tpo2$feature=='pcexon') , 'stop'] <- tpo$exon_chrom_end
  
  # now for introns
  # need to go stop position of previous exon + 1 up to start position of next exon - 1
  
  # how many introns do we have?
  n_introns <- length(which(tpo2$type == 'intron'))
  
  # where we will store start positions
  intron_starts <- rep(NA, n_introns)
  
  # where we will store stop positions
  intron_stops <- rep(NA, n_introns)
  
  # for each intron
  for (i in 1:max(n_introns)) {
    
    # we are filling row number
    rown <- which(tpo2$feature=='intron')[i]
    
    # start position of that intron = stop position of whatever is in the previous row
    intron_starts[i] <- tpo2[rown-1 , 'stop']
    
    # stop position of that intron = start position of whatever is in next row
    intron_stops[i] <- tpo2[rown+1 , 'start']
    
  }
  
  # we place them in tpo2
  tpo2 [which(tpo2$type=='intron') , 'start'] <- intron_starts
  tpo2 [which(tpo2$type=='intron') , 'stop'] <- intron_stops
  
  
  # convert to distances prior to scaling down ------------------------------
  
  cat('\t \t \t \t >>> Convert to intron/exon distances... \n')
  
  tpo2$length <- tpo2$stop - tpo2$start
  

  
  # scale down exons/introns ------------------------------------------------
  
  cat('\t \t \t \t >>> Scale down introns/exons... \n')
  
  # new column to store the scaled down lengths
  tpo2$length_down <- NA
  
  # scale down the exons
  tpo2[which(tpo2$type=='exon'), 'length_down'] <- tpo2[which(tpo2$type=='exon'), 'length'] * scale_exons
  
  # scale down the introns
  tpo2[which(tpo2$type=='intron'), 'length_down'] <- tpo2[which(tpo2$type=='intron'), 'length'] * scale_introns
  
  # rounding will make things more intuitive
  tpo2$length_down <- round(tpo2$length_down)
  
  
  
  # convert back to positions -----------------------------------------------
  
  cat('\t \t \t \t >>> Convert back to positions... \n')
  
  # e.g. length_down are 9, 4, 6
  # start positions should be 0, 9, 13
  # stop positions should be 9, 13, 19
  # >> what we want is cumulative sum
  pos <- c(0, cumsum(tpo2$length_down)) # and add 0 in front
  
  # start positions are 1st pos to n-1th pos
  tpo2$start2 <- pos[1:(length(pos)-1)]
  
  # stop positions are 2nd pos to nth pos
  tpo2$stop2 <- pos[2:length(pos)]
  
  
  
  # extend a little bit the UTRs --------------------------------------------
  
  # better for aesthetics to have the UTRs extend a little bit below the exon
  # ! only the UTRs which have a boundary with a protein-coding exon
  # if exon entirely UTR, should not change its positions
  
  # ! may cause issues if the exon after/before is very small
  # may need to calibrate better in the future
  
  extendby <- 2
  
  ### do it for left UTR ###
  
  # find transition
  # `if` below looks for successive rows where feature transitions leftutr >> pcexon
  # in other words, it looks for left UTR followed by protein-coding exon
  row2change <- c()
  
  for (i in 1:nrow(tpo2)) {
    
    # careful not to go after last feature, can happen if there is no UTR whatsoever
    if(i+1 > nrow(tpo2)) next()
    
    if (tpo2$feature[i]=='leftutr' & tpo2$feature[i+1]=='pcexon') { row2change <- c(row2change, i) }
  }
  
  # there can only be 0 or 1 row2change
  # typically 1, i.e. there is one transition leftutr >> pcexon on the left
  # 0 if no transition, i.e. there is no leftUTR
  if(! length(row2change) %in% c(0, 1)) stop('\t \t \t \t Error: Something unexpected when finding transitions between left UTR and first protein-coding exon \n')
  
  # extend a little bit to the right the last left UTR
  if (length(row2change) > 0) { # if there is a transition left UTR >> protein-coding exon on the left
    tpo2[row2change, 'stop2'] <- tpo2[row2change[1], 'stop2'] + extendby
  }
  
  
  ### do it for right UTR ###
  
  # find transition
  # `if` below looks for successive rows where feature transitions pcexon >> rightutr
  row2change <- c()
  
  for (i in 1:nrow(tpo2)) {
    
    # careful not to go after last feature, can happen if there is no UTR whatsoever
    if(i+1 > nrow(tpo2)) next()
    
    if (tpo2$feature[i]=='pcexon' & tpo2$feature[i+1]=='rightutr') { row2change <- c(row2change, i+1) }
  }
  
  # there can only be 0 or 1 row2change
  # typically 1, i.e. there is one transition leftutr >> pcexon on the left
  # 0 if no transition, i.e. there is no leftUTR
  if(! length(row2change) %in% c(0, 1)) stop('\t \t \t \t Error: Something unexpected when finding transitions between last protein-coding exon and right UTR \n')
  
  # extend a little bit to the left the first right UTR
  if(length(row2change) > 0) { # if there is a transition protein-coding exon >> UTR on the right
    tpo2[row2change, 'start2'] <- tpo2[row2change, 'start2'] - extendby
  }
  

  
  # add target arrows -------------------------------------------------------
  
  # skip if user does not want to mark any position
  
  if (!is.na(positions_tomark[1])) {
    
    # need to find the position to mark and convert its position in new scale
    cat('\t \t \t \t >>> Searching for the position(s) to mark... \n')
    
    # marks will store the positions in new scale
    marks <- c()
    
    for (x in 1:length(positions_tomark)) { # for each position to mark
      
      # clear any mrk created when looking for the previous position to mark
      if(exists('mrk')) {
        rm(mrk)
      }
      
      pampos <- positions_tomark[x]
      
      # try each exon one by one
      
      # how many exons do we need to try:
      npcex <- length(which(tpo2$feature=='pcexon'))
      # ! note, currently cannot mark a position not in a protein-coding exon
      
      for (i in 1:npcex) {
        
        # which row are we looking at?
        rowi <- which(tpo2$feature=='pcexon')[i]
        
        # is position in that exon?
        if (pampos >= tpo2[rowi, 'start'] & pampos <= tpo2[rowi, 'stop']) {
          # >/< or equal as PAM could be exactly the first or last bp of this exon
          
          # if yes, record its position as ratio of total length of that exon
          posr <- which(tpo2[rowi, 'start']:tpo2[rowi, 'stop'] == pampos) / (tpo2[rowi, 'stop'] - tpo2[rowi, 'start'])
          
          # multiply length_down of that exon by position in ratio, will give distance on that exon where the target is
          posdi <- tpo2[rowi, 'length_down'] * posr
          
          # add this distance to start of exon will give where to put the mark on the plot
          mrk <- tpo2[rowi, 'start2'] + posdi
          
        }
        
      }
      
      # check we found the target exon for this position
      
      if(exists('mrk')) {
        cat('\t \t \t \t \t >>> Found site to mark for position #', x, '\n')
      }
      
      if(!exists('mrk')) stop('\n \n \t \t \t \t >>> Error: could not find the target exon for position # ', x, '
                              \n \t \t \t \t Please check in UCSC genome browser or Ensembl exon table that this position falls in an exon.
                              \n \t \t \t \t Note: marking something other than a protein-coding exon (e.g. a UTR or an intron) is not currently supported. \n')
      
      marks <- c(marks, mrk)
      
    }
    
    # turn the marks into a dataframe for ggplot
    mkdf <- as.data.frame(matrix(nrow=length(marks), ncol=2))
    colnames(mkdf) <- c('num', 'plotpos')
    mkdf$num <- 1:nrow(mkdf)
    mkdf$plotpos <- marks
    
  }

  
  
  
  # ready to plot -----------------------------------------------------------
  
  cat('\t \t \t \t >>> Plot... \n')
  
  # we are plotting intron as simply one line, and exons plotted a second layer on top
  # so only take the exon positions for plotting
  tpo2ex <- subset(tpo2, type=='exon')
  
  map <- ggplot(tpo2ex) +
    geom_rect(aes(xmin=min(start2), xmax=max(stop2), ymin=(10-intron_height/2), ymax=(10+intron_height/2)), fill=intron_colour) +
    
    { if(!roundCorners) geom_rect(aes(xmin=start2, xmax=stop2, ymin=ystart, ymax=ystop), fill=exon_colour) } +
    
    { if(roundCorners) statebins:::geom_rrect(aes(xmin=start2, xmax=stop2, ymin=ystart, ymax=ystop),
                                              radius=grid::unit(corner_radius, 'pt'), fill=exon_colour) } +
    
    # add arrows for marked positions
    # arrowsize / 15 is by trial and error here
    # need to leave enough space for arrow to grow if user puts a bigger arrowsize, but at the same time small enough to hide the line below the arrowhead
    { if (!is.na(positions_tomark[1]) & arrowpos=='above') # skip if user does not want to mark any position
      geom_segment(data=mkdf, aes(x=plotpos, y=10+(exon_height/2)+(arrowsize/15), xend=plotpos, yend=10+(exon_height/2)),
                 lineend='butt',
                 linejoin='mitre',
                 arrow = arrow(length=unit(arrowsize, 'mm'),
                               type='closed'),
                 colour='#EE7163') } +
    
    { if (!is.na(positions_tomark[1]) & arrowpos=='below') # skip if user does not want to mark any position
      geom_segment(data=mkdf, aes(x=plotpos, y=10-(exon_height/2)-(arrowsize/15), xend=plotpos, yend=10-(exon_height/2)),
                   lineend='butt',
                   linejoin='mitre',
                   arrow = arrow(length=unit(arrowsize, 'mm'),
                                 type='closed'),
                   colour='#EE7163') } +
    
    theme_minimal() +
    theme(
      panel.grid=element_blank(),
      legend.position='none',
      axis.text=element_blank(),
      axis.title=element_blank(),
      plot.margin=margin(t=0, r=-6, b=0, l=-11)
    )
  
  # export plot -------------------------------------------------------------
  
  cat('\t \t \t \t >>> Export plot as', paste0(exportFolder, '/', transcript_id, '.pdf'), '... \n')
  
  ggsave(paste0(exportFolder, '/', transcript_id, '.pdf'), plot=map, width=width, height=height, units='mm')
  
  

  # return plot -------------------------------------------------------------
  # so it displays
  return(map)
  
  
}

