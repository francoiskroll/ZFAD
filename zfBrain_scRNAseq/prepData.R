# exploring single-cell RNAseq data from Raj et al., 2020 (Schier lab)
# here: downloading/preparing the data
# beware -- takes a while to run and you need ~ 90 Gb free on your disk

### note 27/04/2023: updated paths but did not run again below


# packages ----------------------------------------------------------------

library(here)
library(RCurl)
library(openxlsx)
library(GEOquery)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)


# download cluster info ---------------------------------------------------

# before actual data; cluster info is in ZFBrainAtlasMaster.tsv
# here https://github.com/brlauuu/zf_brain/blob/master/data/ZFBrainAtlasMaster.tsv
cls <- read.csv(text=getURL('https://raw.githubusercontent.com/brlauuu/zf_brain/master/data/ZFBrainAtlasMaster.tsv'),
                sep='\t') # clusters
# (I use read.csv but it is not actually a csv but a tsv, hence sep='\t')

# it imports a column X which is just row number, delete it
cls[,1] <- NULL

# fix column names
colnames(cls) <- c('stage', 'cluster', 'match_previousstage', 'celltype', 'enrichedmarkers')

# will also save it locally
write.csv(cls, here('zfBrain_scRNAseq', 'zfclusters.csv'), row.names=FALSE)

# first few stages do not have the same name in cluster information and file name below
# correspondence (from https://zfin.org/zf_info/zfbook/stages/index.html)
# 6 somites = 12 hpf
# 10 somites = 14 hpf
# 14 somites = 16 hpf
# 18 somites = 18 hpf


# download Seurat data ----------------------------------------------------

dir.create(here('zfBrain_scRNAseq', 'seuratdata'))

# I wrote info and download urls of files in filesMeta.xlsx
# from Raj 2020 GEO repository GSE158142
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158142
# to get http/ftp url:
  # right click http or ftp link > Copy Link Address
filesmeta <- read.xlsx(here('zfBrain_scRNAseq', 'filesMeta.xlsx'))

# will store Seurat object in a list
# one Seurat object per developmental stage, 12 developmental stages in total
seus <- vector(length=12, mode='list') # Seurat objects

# call each element by age
names(seus) <- filesmeta$age[1:12]

# ! by default download_file() will timeout after 60sec, but we need longer to download the bigger files
getOption('timeout')
options(timeout=60*60*24) # 24hours

# for each developmental stage:
# download the corresponding Seurat object locally
# unzip it
# import the RDS file into R
# add it to the list

for(i in filesmeta$age_order) { # i goes from 1 to 12
  # take the link from the meta info file
  urlSeu <- filesmeta[which(filesmeta$age_order==i), 'fttp_url']
  fileSeu <- filesmeta[which(filesmeta$age_order==i), 'file']
  
  # download the file
  download.file(url=urlSeu,
                destfile=here('zfBrain_scRNAseq', 'seuratdata', fileSeu))
  
  # unzip it
  gunzip(here('zfBrain_scRNAseq', 'seuratdata', fileSeu))
  
  # ! file name will not have .gz anymore; so remove last 3 characters
  fileSeu <- substr(fileSeu, 1, nchar(fileSeu)-3)
  
  # import the RDS into R + put the object in the list
  seus[[i]] <- UpdateSeuratObject(readRDS(here('zfBrain_scRNAseq', 'seuratdata', fileSeu)))
  # UpdateSeuratObject() looks important as in global.R here https://github.com/brlauuu/zf_brain
}

# takes a while to download everything, so save the list for later
saveRDS(object=seus, file=here('zfBrain_scRNAseq', 'seuratAll.rds'))

# could delete folder seuratdata now, everything is stored in seuratAll.rds

# we are mostly interested in 5-dpf data, so save that dataset only as a separate file
# and we can put seuratAll.rds on Dropbox "online-only" most of the time
seu5dpf <- seus$dpf5
rm(seus)
gc()
saveRDS(object=seu5dpf, file=here('zfBrain_scRNAseq', 'seurat5dpf.rds'))
