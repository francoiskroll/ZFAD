#####################################################
# ~ FbyFcomparison: function getSleepFromMat ~
#
# function to get sleep parameters from mat file
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################



# packages ----------------------------------------------------------------

library(R.matlab)
library(here)
library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
library(naturalsort)



# function getSleepFromMat ------------------------------------------------

# parameter: sleep parameter to extract
# options are identical to FramebyFrame parameters:
# sleepHours, sleepNumNaps, sleepNapDuration, sleepLatency

getSleepFromMat <- function(matPath,
                            parameter,
                            boxnum=1) {

  ### import mat file
  mat <- readMat(matPath)$geno

  # how many groups are there?
  ngrps <- ncol(mat[,,1]$name)
  cat('\t \t \t \t >>> Found', ngrps,'groups in .mat file \n')

  ### get fish IDs
  fids <- unlist(sapply(1:ngrps, function(ng) {
    unlist(mat[,,1]$fishID[,ng])
  }))
  # if e.g. there are two groups
  # this is simply all larvae (well numbers) in group1, then all larvae (well numbers) in group2
  # simply in the order given in genotype file

  # ! if boxnum=2, need to add 96
  if(boxnum==2) {
    fids <- fids+96
  }

  # add f in front of IDs so easier to match with FramebyFrame parameter table
  fids <- sprintf('f%i', fids)

  # import each group's data in a list
  # each slot of the list is datapoints for all larvae of one group


  ### get parameter

  ##########
  if(parameter=='sleepHours') {

    ## import DAY data
    palday <- lapply(1:ngrps, function(ng) {
      as.data.frame( t(as.data.frame(mat[,,1]$summarytable[,,1]$sleep[,,1]$day[[ng]])) ) # an ugly line...
    })
    # append the list, so we have all larvae of group1, then all larvae of group2 below, etc.
    paday <- as.data.frame(rbindlist(palday))
    colnames(paday) <- sprintf('day%i', 0:(ncol(paday)-1)) # so day0, day1, ...

    ## import NIGHT data
    palnight <- lapply(1:ngrps, function(ng) {
      as.data.frame( t(as.data.frame(mat[,,1]$summarytable[,,1]$sleep[,,1]$night[[ng]])) )
    })
    # append the list, so we have all larvae of group1, then all larvae of group2 below, etc.
    panight <- as.data.frame(rbindlist(palnight))
    colnames(panight) <- sprintf('night%i', 0:(ncol(panight)-1)) # so night0, night1, ...

    ## merge day & night data
    pa <- cbind(paday, panight)

    # total sleep currently in minutes, we want in hours
    pa <- pa/60

    ##########
  } else if(parameter=='sleepNumNaps') {

    ## import DAY data
    palday <- lapply(1:ngrps, function(ng) {
      as.data.frame( t(as.data.frame(mat[,,1]$summarytable[,,1]$sleepBout[,,1]$day[[ng]])) ) # an ugly line...
    })
    # append the list, so we have all larvae of group1, then all larvae of group2 below, etc.
    paday <- as.data.frame(rbindlist(palday))
    colnames(paday) <- sprintf('day%i', 0:(ncol(paday)-1)) # so day0, day1, ...

    ## import NIGHT data
    palnight <- lapply(1:ngrps, function(ng) {
      as.data.frame( t(as.data.frame(mat[,,1]$summarytable[,,1]$sleepBout[,,1]$night[[ng]])) )
    })
    # append the list, so we have all larvae of group1, then all larvae of group2 below, etc.
    panight <- as.data.frame(rbindlist(palnight))
    colnames(panight) <- sprintf('night%i', 0:(ncol(panight)-1)) # so night0, night1, ...

    ## merge day & night data
    pa <- cbind(paday, panight)

    ##########
  } else if(parameter=='sleepNapDuration') {

    ## import DAY data
    palday <- lapply(1:ngrps, function(ng) {
      as.data.frame( t(as.data.frame(mat[,,1]$summarytable[,,1]$sleepLengthmean[,,1]$day[[ng]])) ) # an ugly line...
    })
    # append the list, so we have all larvae of group1, then all larvae of group2 below, etc.
    paday <- as.data.frame(rbindlist(palday))
    colnames(paday) <- sprintf('day%i', 0:(ncol(paday)-1)) # so day0, day1, ...

    ## import NIGHT data
    palnight <- lapply(1:ngrps, function(ng) {
      as.data.frame( t(as.data.frame(mat[,,1]$summarytable[,,1]$sleepLengthmean[,,1]$night[[ng]])) )
    })
    # append the list, so we have all larvae of group1, then all larvae of group2 below, etc.
    panight <- as.data.frame(rbindlist(palnight))
    colnames(panight) <- sprintf('night%i', 0:(ncol(panight)-1)) # so night0, night1, ...

    ## merge day & night data
    pa <- cbind(paday, panight)

  } else if(parameter=='sleepLatency') {

    ## import DAY data
    palday <- lapply(1:ngrps, function(ng) {
      as.data.frame( t(as.data.frame(mat[,,1]$summarytable[,,1]$sleepLatency[,,1]$day[[ng]])) ) # an ugly line...
    })
    # append the list, so we have all larvae of group1, then all larvae of group2 below, etc.
    paday <- as.data.frame(rbindlist(palday))
    colnames(paday) <- sprintf('day%i', 0:(ncol(paday)-1)) # so day0, day1, ...

    ## import NIGHT data
    palnight <- lapply(1:ngrps, function(ng) {
      as.data.frame( t(as.data.frame(mat[,,1]$summarytable[,,1]$sleepLatency[,,1]$night[[ng]])) )
    })
    # append the list, so we have all larvae of group1, then all larvae of group2 below, etc.
    panight <- as.data.frame(rbindlist(palnight))
    colnames(panight) <- sprintf('night%i', 0:(ncol(panight)-1)) # so night0, night1, ...

    ## merge day & night data
    pa <- cbind(paday, panight)

  }


  ### add fish IDs as a column of parameter table
  # and parameter we chose
  if(length(fids)!=nrow(pa)) stop('\t \t \t \t >>> Error: not the same number of fish IDs than rows in parameter table \n')

  pa <- pa %>%
    add_column(fish=fids, .before=1) %>%
    add_column(parameter=parameter, .before=1)


  ### order by increasing fish ID
  # ! not in anymore simply all larvae grp1 then all larvae grp2 below
  pa <- pa[naturalorder(pa$fish),]


  ### return pa
  return(pa)

}
