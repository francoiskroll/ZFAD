# write Forward_Reads and Reverse_Reads

library(openxlsx)



fqnms <- list.files(path='~/Dropbox/phd/220524_miseq/data/reads/filterfastq/')

# get gene names
genms <- sapply(strsplit(fqnms, split='_'), function(nm) { nm[2] } )

# order of genes you want
genorder <- c('apoea', 'apoeb', 'appa', 'appb', 'sorl1', 'cd2ap', 'clu')

# order fastq filenames in same order as samples in plate
fqnms <- fqnms[order(match(genms, genorder))]

# preallocate dataframe
fqdf <- as.data.frame(matrix(nrow=length(fqnms)/2, ncol=4))
colnames(fqdf) <- c('Forward_Reads', 'Reverse_Reads', 'Group', 'Control')

# put all the R1 files in first column, all the R2 files in second column
fqdf$Forward_Reads <- fqnms[seq(1, length(fqnms), 2)] # all even indices, should be all R1
fqdf$Reverse_Reads <- fqnms[seq(2, length(fqnms), 2)] # all odd indices, should be all R2

fqdf$Group <-  rep(1:17, each=8)

fqdf$Control <- rep(c(1, 1, 0, 0, 0, 0, 0, 0), 17)

# write
write.xlsx(x=fqdf,
           file='~/Dropbox/phd/220524_miseq/data/readnames.xlsx',
           overwrite=TRUE)
