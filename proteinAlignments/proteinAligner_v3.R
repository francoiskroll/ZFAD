# proteinAligner
# draws an aligment schematic to compare human protein with zebrafish protein

# essentially drawing a schematic an alignment in https://www.ebi.ac.uk/Tools/msa/clustalo/
# colours correspond to the symbols . : * in the alignment

# find full example of function call at the end


### version history ###

# v2:

# v1: if a similarity group is not present (e.g. there is 0 "highly similar" match), the colours were getting shifted
# now corrected


### ###

# packages & small functions ----------------------------------------------

# BiocManager::install("msa")

library(msa)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)

# function matchSimGroups(...)
# small function to find for one amino acid all its matches in groups of strongly similar amino acids (ssim below) or weakly similar amino acids (wsim below)

# aa = a given amino acid, e.g. A
# whichsim = list ssim or wsim below
matchSimGroups <- function(aa, whichsim) {
  
  # where we will store the matches
  grpmats <- c() # group matches (in positions)
  
  sapply(1:length(whichsim), function(gp) { # gp for group position (i.e. position in the list)
    if (aa %in% whichsim[[gp]]) {grpmats <<- c(grpmats, gp)}
  })
  
  # return the matches for that amino acid
  return(grpmats)
  
}


# similarity groups -------------------------------------------------------
# source: https://www.ebi.ac.uk/seqdb/confluence/display/THD/Help+-+Clustal+Omega+FAQ#HelpClustalOmegaFAQ-Whatdotheconsensussymbolsmeaninthealignment?

# will store groups in a list of vectors
# function matchSimGroups(...) will then scroll through the groups in a list and return the ones that include the amino acid we are looking for

# groups of amino acids with strongly similar properties
ssim <- vector(mode='list', length=9) # ssim for strongly similar
ssim[[1]] <- c('S', 'T', 'A')
ssim[[2]] <- c('N', 'E', 'Q', 'K')
ssim[[3]] <- c('N', 'H', 'Q', 'K')
ssim[[4]] <- c('N', 'D', 'E', 'Q')
ssim[[5]] <- c('Q', 'H', 'R', 'K')
ssim[[6]] <- c('M', 'I', 'L', 'V')
ssim[[7]] <- c('M', 'I', 'L', 'F')
ssim[[8]] <- c('H', 'Y')
ssim[[9]] <- c('F', 'Y', 'W')

# groups of amino acids with weakly similar properties
wsim <- vector(mode='list', length=9) # wsim for weakly similar
wsim[[1]] <- c('C', 'S', 'A')
wsim[[2]] <- c('A', 'T', 'V')
wsim[[3]] <- c('S', 'A', 'G')
wsim[[4]] <- c('S', 'T', 'N', 'K')
wsim[[5]] <- c('S', 'T', 'P', 'A')
wsim[[6]] <- c('S', 'G', 'N', 'D')
wsim[[7]] <- c('S', 'N', 'D', 'E', 'Q', 'K')
wsim[[8]] <- c('N', 'D', 'E', 'Q', 'H', 'K')
wsim[[8]] <- c('N', 'E', 'Q', 'H', 'R', 'K')
wsim[[9]] <- c('F', 'V', 'L', 'I', 'M')
wsim[[10]] <- c('H', 'F', 'Y')

# one amino acid is often included in multiple groups
# so we will need to find the groups which include amino acid 1; 
# for this, we will need to check if intersect is more than 1, i.e. are the amino acids present together in at least one group
# e.g.
# length(intersect(matchSimGroups(aa='K', whichsim=wsim), matchSimGroups(aa='S', whichsim=wsim)))


# function proteinAligner(...) --------------------------------------------

# human = human protein sequence (reference to align the zebrafish sequence to), e.g. 'TLALLAIFKKALPALPI'
# zebrafish = zebrafish protein sequence, e.g. 'LLLLAIFKALPALP'

# cols = colours for alignment scores; order is 0, 1, 2, 3, 4, 5
# meaning is
# 0 = gap, should typically be white or transparent
# 1 = unknown (i.e. in front of gap)
# 2 = mismatch
# 3 = weakly similar
# 4 = highly similar
# 5 = match

# 6 = human (match)

# about defaults:
# mismatch is red = #cb2a20
# match is grey blue = #697a87
# need two intermediary colours for weakly similar (closer to mismatch so more red) and highly similar (closer to match so more greyblue)
# to get back I used https://www.colorhexa.com/cb2a20-to-697a87

# proteinname, e.g. 'psen1'
# only used for file name of pdf plot

# height = height of pdf plot in mm
# width = width of pdf plot in mm

# exportFolder = where to save the plot

proteinAligner <- function(human,
                           zebrafish,
                           labelhumanpos=NA,
                           cols=c('#ffffff', '#97a4ae', '#cb2a20', '#a24b4b', '#8a5f65', '#697a87', '#595E60'), 
                           proteinname,
                           spaceBetween=0,
                           height=90,
                           width=30,
                           exportFolder) {
  
  if (length(cols)!=7) stop('\t \t \t \t >>> Error, not 7 colours for alignment scores \n')
  
  
  # remove any return/new line character in the sequences -------------------
  
  human <- gsub("[\r\n]", '', human)
  zebrafish <- gsub("[\r\n]", '', zebrafish)
  
  
  # align the two sequences -------------------------------------------------
  
  ali <- msa(c(human, zebrafish), type='protein', method='ClustalOmega') # ali for alignment
  # method can be changed, here ClustalOmega
  
  print(ali)
  
  # put the aligned sequences in a dataframe
  alm <- as.data.frame(matrix(nrow=2, ncol=nchar(as.character(unmasked(ali)[1]))))
  # alm for alignment matrix, simply the alignment in a more convenient table we can loop through
  
  # place the human sequence in the alignment matrix
  alm[1,] <- strsplit(as.character(unmasked(ali)[1]), split='')[[1]]
  
  # place the zebrafish sequence in the alignment matrix
  alm[2,] <- strsplit(as.character(unmasked(ali)[2]), split='')[[1]]
  
  
  # preallocate protein comparison matrix -----------------------------------
  
  # matrix will store 'alignment score' for each position
  # options are, broadly from worse alignment to best alignment;
  # 0 = gap
  # 1 = unknown (i.e. in front of gap)
  # 2 = mismatch
  # 3 = weakly similar
  # 4 = highly similar
  # 5 = match
  
  pcm <- as.data.frame(matrix(nrow=nchar(as.character(unmasked(ali)[1])), ncol=3)) # protein comparison matrix
  colnames(pcm) <- c('pos', 'human', 'zebrafish')
  pcm$pos <- 1:nrow(pcm)
  
  # add original positions in human sequence
  # this is simply 1, 2, 3, ..., total length of human sequence
  # but should not count gaps, set to 0
  # we can tell from alm where the gaps are (if any)
  
  # preallocate human positions after alignment
  hpos <- c()
  
  # initiate position counter at 0
  # so that first position is 1
  lastpos <- 0
  
  # then fill in position by position
  for(i in 1:length(alm[1,])) {
    
    # if position is gap, we add 0
    if(alm[1,i] == '-') {
      hpos <- c(hpos, 0)
      
    # if not, we add last position + 1
    } else {
      hpos <- c(hpos, lastpos + 1)
      # and we keep track in lastpos
      lastpos <- lastpos + 1
    }
    
  }
  
  # maximum human position should the total length of the human sequence
  if(max(hpos) != nchar(human)) stop('\t \t \t \t >>> Error when annotating original human positions in aligment. \n')
  
  # add human positions to pcm
  pcm <- pcm %>%
    mutate(hpos=hpos, .after='pos')
  
  # human sequence is reference
  # so it is never incorrect
  # only possibilities are match (score 6, which is only for human) or gap (score 0)
  # can put 6 for every position, then replace by 0 where needed below
  pcm$human <- 6
  
  
  # fill in protein comparison matrix ---------------------------------------
  
  # note about below; a given `else if` is only triggered if previous `else if` (or `if`) was not triggered
  # so should be read as: "try full match, if does not work try gap in human sequence, if does not work try gap in zebrafish sequence", etc.
  # hence the order of similarity (score 5, 4, 3, 2), as it should be trying full match first, strongly similar second, weakly similar third, mismatch fourth
  # or in other words, last resort (last `else`) is complete mismatch
  
  sapply(1:ncol(alm), function(pos) {
    
    pa <- alm[,pos] # alignment at this position
    
    # then list all possibilities
    # only changing zebrafish sequence, except when gap in human sequence
    
    # if gap in human sequence
    if (pa[1] == '-') {
      pcm[pos, 'human'] <<- 0 # gap in human sequence
      pcm[pos, 'zebrafish'] <<- 1 # unknown in zebrafish sequence
    }
    
    # if gap in zebrafish sequence
    else if (pa[2] == '-') {
      pcm[pos, 'zebrafish'] <<- 0 # gap in zebrafish sequence
    }
    
    # I think gap in both cannot exist
    
    # if full match
    else if (pa[1] == pa[2]) {
      pcm[pos, 'zebrafish'] <<- 5 # match in zebrafish sequence
    }
    
    # if strongly similar
    # i.e. there is more than one group of strongly similar amino acids where both the zebrafish amino acid and the human amino acid are present
    # matchSimGroups for human amino acid finds which similarity group(s) (if any) includes it
    # matchSimGroups for zebrafish amino acid finds which similarity group(s) (if any) includes it
    # intersect between the two gives any groups that include both amino acids
    # length > 1 is just to check if there is at least one group that include both amino acids
    else if ( length(intersect(matchSimGroups(aa=pa[1], whichsim=ssim), matchSimGroups(aa=pa[2], whichsim=ssim))) > 0 ) {
      pcm[pos, 'zebrafish'] <<- 4 # strongly similar in zebrafish sequence
    }
    
    # if weakly similar
    # i.e. there is more than one group of weakly similar amino acids where both the zebrafish amino acid and the human amino acid are present
    # logic same as above
    else if ( length(intersect(matchSimGroups(aa=pa[1], whichsim=wsim), matchSimGroups(aa=pa[2], whichsim=wsim))) > 0 ) {
      pcm[pos, 'zebrafish'] <<- 3 # weakly similar in zebrafish sequence
    }
    
    # if mismatch (not similar)
    # i.e. if none of the options above were triggered, last resort is mismatch
    else {
      pcm[pos, 'zebrafish'] <<- 2 # mismatch in zebrafish sequence
    }
    
  })
  
  

# return identity % -------------------------------------------------------
  
  # it is simply number of 5s divided by total length
  cat('\t \t \t \t >>> Identity:', length(which(pcm$zebrafish==5)), '/', nrow(pcm), '=',
      round(100 * (length(which(pcm$zebrafish==5)) / nrow(pcm)), 2), '% \n')
  
  

# return similarity % -----------------------------------------------------
  # disclaimer: I have not read much about how it is typically defined
  # will simply define here as equal or highly similar
  cat('\t \t \t \t >>> Similarity (identical or highly similar):', length(which(pcm$zebrafish %in% c(4,5))), '/',
      nrow(pcm), '=',
      round(100 * (length(which(pcm$zebrafish %in% c(4, 5))) / nrow(pcm)), 2), ' % \n')
  
  
  ### if need to label some human positions
  if(!is.na(labelhumanpos[1])) {
    
    # we simply swap those positions to mismatch in pcm
    pcm[which(pcm$hpos %in% labelhumanpos), 'human'] <- 2
    
  }

# plot --------------------------------------------------------------------
  
  # add names to the colour vector so we make sure correct colour is used for correct similarity level
  # corrects v1 issue: when a similarity level was never present (e.g. no "highly similar"), it was not matching the colours correctly
  names(cols) <- 0:6
  
  pcm$zebrafish <- factor(pcm$zebrafish, levels=c(0:6))
  pcm$human <- factor(pcm$human, levels=c(0:6)) # human can only be match (6) or gap (0)
  
  # protein comparison plot
  pcp <- ggplot(pcm) +
    geom_rect(aes(xmin=pos, xmax=pos+1, ymin=13, ymax=14, fill=human)) +
    geom_rect(aes(xmin=pos, xmax=pos+1, ymin=8-spaceBetween, ymax=12-spaceBetween, fill=zebrafish)) +
    scale_fill_manual(values=cols) +
    theme_minimal() +
    theme(
      axis.title=element_blank(),
      axis.text=element_blank(),
      panel.grid=element_blank(),
      legend.position='none',
      plot.margin=margin(t=-2, r=-8, b=-3, l=-14)
    )
  
  
# save plot ---------------------------------------------------------------

  
  cat('\n \n \t \t \t \t >>> Export plot as', paste0(exportFolder, '/', proteinname, '.pdf'), '... \n')
  
  ggsave(paste0(exportFolder, '/', proteinname, '.pdf'), plot=pcp, width=width, height=height, units='mm')
  
  

# return plot -------------------------------------------------------------
  # so it displays in RStudio
  
  return(pcp)
  
}


# full example ------------------------------------------------------------

# can be copy-pasted directly in another script

# source('~/Dropbox/phd/utilities/proteinAligner.R')
#
# # from https://www.uniprot.org/uniprot/P49768.fasta
# human_psen1 <- 'MTELPAPLSYFQNAQMSEDNHLSNTVRSQNDNRERQEHNDRRSLGHPEPLSNGRPQGNSR
# QVVEQDEEEDEELTLKYGAKHVIMLFVPVTLCMVVVVATIKSVSFYTRKDGQLIYTPFTE
# DTETVGQRALHSILNAAIMISVIVVMTILLVVLYKYRCYKVIHAWLIISSLLLLFFFSFI
# YLGEVFKTYNVAVDYITVALLIWNFGVVGMISIHWKGPLRLQQAYLIMISALMALVFIKY
# LPEWTAWLILAVISVYDLVAVLCPKGPLRMLVETAQERNETLFPALIYSSTMVWLVNMAE
# GDPEAQRRVSKNSKYNAESTERESQDTVAENDDGGFSEEWEAQRDSHLGPHRSTPESRAA
# VQELSSSILAGEDPEERGVKLGLGDFIFYSVLVGKASATASGDWNTTIACFVAILIGLCL
# TLLLLAIFKKALPALPISITFGLVFYFATDYLVQPFMDQLAFHQFYI'
# 
# # from https://www.uniprot.org/uniprot/Q9W6T7.fasta
# zebrafish_psen1 <- 'MADLVQNAANNVLNDGMDTSRHTSSTAAPPSRNEVELNGQPPTAPPPQVVTDSEEDEDEE
# LTLKYGAKHVIMLFIPVTLCMVVVVATIKSVSFYTQKDGQQLIYTPFREDTETVGQRALH
# SMLNAIIMISVIVVMTLVLVVLYKYRCYKVIQAWLFFSNLLLLFFFSLIYLGEVFKTYNV
# AMDYFTLALIIWNFGVVGMICIHWKGPLRLQQAYLIMISALMALVFIKYLPEWTAWLILA
# AISVYDLLAVLCPKGPLRILVETAQERNEAIFPALIYSSTMVWLFNMADSAETRNNSSHP
# VPQQENQVVAMAPTAQAEDDGGFTPAWVDHQQHQLGPMQSTEESRRQIQEMPSARPPPPA
# DDDEERGVKLGLGDFIFYSMLVGKASATASGDWNTTLACFVAILIGLCLTLLLLAIFKKA
# LPALPISITFGLVFYFATDNLVRPFMDQLAVHQFYI'
# 
# proteinAligner(human=human_psen1,
#                zebrafish=zebrafish_psen1,
#                cols=c('#ffffff', '#97a4ae', '#cb2a20', '#a24b4b', '#8a5f65', '#697a87'),
#                proteinname='psen1',
#                width=90,
#                height=30,
#                exportFolder='~/Dropbox')
