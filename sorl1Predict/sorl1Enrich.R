###################################################
# ~~~ ZFAD: predictive pharmacology ~~~

# test with sorl1 fingerprint

# Francois Kroll 2022
# francois@kroll.be
###################################################


# packages/functions ------------------------------------------------------

library(here)
library(tictoc)

source(here('annotateDrugDb', 'drawEnrich_v4.R'))
source(here('annotateDrugDb', 'ggEnrich.R'))
source(here('annotateDrugDb', 'legacyFingerprint', 'legacyFingerprint.R'))
source(here('annotateDrugDb', 'legacyFingerprint', 'gglegacyFingerprint.R'))


# calculate fingerprint ---------------------------------------------------

# calculate legacy fingerprint
lfp <- legacyFingerprint(matPath=here('220531_SORL1', 'legacyMiddur', '220531_14', '220531_14.mat'),
                         conGrp='scr',
                         treGrp='sorl1',
                         days=c(2,3),
                         nights=c(2,3))


# rank drugDb -------------------------------------------------------------

# rank drug db in comparison to fingerprint
vdbr <- rankDrugDb(legacyFgp=lfp,
                   dbPath=here('annotateDrugDb', 'drugDb.csv'), 
                   metric='cosine')
write.csv(vdbr, here('predPharma', 'sorl1Enrich_f0exp1', 'sorl1Ranked_f0exp1.csv'), row.names=FALSE)


# enrichment of indications -----------------------------------------------

tic()
inr <- drugEnrichment(vdbr=vdbr,
                      namesPath=here('annotateDrugDb', 'compounds.csv'),
                      annotationPath=here('annotateDrugDb', 'TTDindications.csv'),
                      annotation='indications',
                      whichRank='rankeq',
                      minNex=3,
                      ndraws=100000,
                      alphaThr=0.05,
                      statsExport=here('predPharma', 'sorl1Enrich_f0exp1', 'sorl1Indications_f0exp1.csv'))
toc()

# 3 remain significant after BenHoch correction
# Depression (top1)
ggDraws(vdbr=vdbr,
        namesPath=here('annotateDrugDb', 'compounds.csv'),
        annotationPath=here('annotateDrugDb', 'TTDindications.csv'),
        annotation='indications',
        testAnnotation='Depression',
        ndraws=100000,
        whichRank='rankeq',
        exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'depression.pdf'),
        width=100,
        height=100)

ggBarcode(vdbr=vdbr,
          namesPath=here('annotateDrugDb', 'compounds.csv'),
          annotationPath=here('annotateDrugDb', 'TTDindications.csv'),
          annotation='indications',
          testAnnotation='Depression',
          barwidth1=2,
          barwidth2=25,
          exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'bc_depression.pdf'),
          width=200,
          height=100)

# Skin imperfections (top2)
ggDraws(vdbr=vdbr,
        namesPath=here('annotateDrugDb', 'compounds.csv'),
        annotationPath=here('annotateDrugDb', 'TTDindications.csv'),
        annotation='indications',
        testAnnotation='Skin imperfections',
        ndraws=100000,
        whichRank='rankeq',
        exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'skinimperfections.pdf'),
        width=100,
        height=100)

ggBarcode(vdbr=vdbr,
          namesPath=here('annotateDrugDb', 'compounds.csv'),
          annotationPath=here('annotateDrugDb', 'TTDindications.csv'),
          annotation='indications',
          testAnnotation='Skin imperfections',
          barwidth1=2,
          barwidth2=25,
          exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'bc_skinimperfections.pdf'),
          width=200,
          height=100)
# Schizophrenia (top3)
ggDraws(vdbr=vdbr,
        namesPath=here('annotateDrugDb', 'compounds.csv'),
        annotationPath=here('annotateDrugDb', 'TTDindications.csv'),
        annotation='indications',
        testAnnotation='Schizophrenia',
        ndraws=100000,
        whichRank='rankeq',
        exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'schizophrenia.pdf'),
        width=100,
        height=100)

ggBarcode(vdbr=vdbr,
          namesPath=here('annotateDrugDb', 'compounds.csv'),
          annotationPath=here('annotateDrugDb', 'TTDindications.csv'),
          annotation='indications',
          testAnnotation='Schizophrenia',
          barwidth1=2,
          barwidth2=25,
          exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'bc_schizophrenia.pdf'),
          width=200,
          height=100)


# enrichment of TTD targets -----------------------------------------------

tic()
ttd <- drugEnrichment(vdbr=vdbr,
                      namesPath=here('annotateDrugDb', 'compounds.csv'),
                      annotationPath=here('annotateDrugDb', 'TTDtargets.csv'),
                      annotation='TTDtargets',
                      whichRank='rankeq',
                      minNex=3,
                      ndraws=100000,
                      alphaThr=0.05,
                      statsExport=here('predPharma', 'sorl1Enrich_f0exp1', 'sorl1TTDtargets_f0exp1.csv'))
toc()

# 8 are significant after BenHoch correction
# top1 is Serotonin transporter (SERT) - SLC6A4
ggDraws(vdbr=vdbr,
        namesPath=here('annotateDrugDb', 'compounds.csv'),
        annotationPath=here('annotateDrugDb', 'TTDtargets.csv'),
        annotation='TTDtargets',
        testAnnotation='T27812',
        ndraws=100000,
        whichRank='rankeq',
        exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'slc6a4.pdf'),
        width=100,
        height=100)

ggBarcode(vdbr=vdbr,
          namesPath=here('annotateDrugDb', 'compounds.csv'),
          annotationPath=here('annotateDrugDb', 'TTDtargets.csv'),
          annotation='TTDtargets',
          testAnnotation='T27812',
          barwidth1=2,
          barwidth2=25,
          exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'bc_slc6a4.pdf'),
          width=200,
          height=100)


# enrichment of KEGG pathways ---------------------------------------------

tic()
keg <- drugEnrichment(vdbr=vdbr,
                      namesPath=here('annotateDrugDb', 'compounds.csv'),
                      annotationPath=here('annotateDrugDb', 'TTDkegg.csv'),
                      annotation='KEGG',
                      whichRank='rankeq',
                      minNex=3,
                      ndraws=100000,
                      alphaThr=0.05,
                      statsExport=here('predPharma', 'sorl1Enrich_f0exp1', 'sorl1KEGG_f0exp1.csv'))
toc()

# 8 are significant after BenHoch correction
# top1 is Serotonergic synapse
ggDraws(vdbr=vdbr,
        namesPath=here('annotateDrugDb', 'compounds.csv'),
        annotationPath=here('annotateDrugDb', 'TTDkegg.csv'),
        annotation='KEGG',
        testAnnotation='hsa04726',
        ndraws=100000,
        whichRank='rankeq',
        exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'serotonergicSynapse.pdf'),
        width=100,
        height=100)

ggBarcode(vdbr=vdbr,
          namesPath=here('annotateDrugDb', 'compounds.csv'),
          annotationPath=here('annotateDrugDb', 'TTDkegg.csv'),
          annotation='KEGG',
          testAnnotation='hsa04726',
          barwidth1=2,
          barwidth2=25,
          exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'bc_serotonergicSynapse.pdf'),
          width=200,
          height=100)

# ironically, one of the worst (pval = 1.0) is Alzheimer's disease
ggDraws(vdbr=vdbr,
        namesPath=here('annotateDrugDb', 'compounds.csv'),
        annotationPath=here('annotateDrugDb', 'TTDkegg.csv'),
        annotation='KEGG',
        testAnnotation='hsa05010',
        ndraws=100000,
        whichRank='rankeq',
        exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'alzheimerDisease.pdf'),
        width=100,
        height=100)

ggBarcode(vdbr=vdbr,
          namesPath=here('annotateDrugDb', 'compounds.csv'),
          annotationPath=here('annotateDrugDb', 'TTDkegg.csv'),
          annotation='KEGG',
          testAnnotation='hsa05010',
          barwidth1=2,
          barwidth2=25,
          exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'bc_alzheimerDisease.pdf'),
          width=200,
          height=100)


# zebrafish STITCH targets ------------------------------------------------

tic()
zstitch <- drugEnrichment(vdbr=vdbr,
                          namesPath=here('annotateDrugDb', 'compounds.csv'),
                          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
                          annotation='zebrafishSTITCH',
                          whichRank='rankeq',
                          minScore=900,
                          minNex=10,
                          ndraws=100000,
                          alphaThr=0.05,
                          maxPval=0.2,
                          statsExport=here('predPharma', 'sorl1Enrich_f0exp1', 'sorl1zSTITCH_f0exp1.csv'))
toc()

# top33 are significant after BenHoch correction
# top annotated one is solute carrier family 6 member 17
ggDraws(vdbr=vdbr,
        namesPath=here('annotateDrugDb', 'compounds.csv'),
        annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
        annotation='zebrafishSTITCH',
        testAnnotation='ENSDARP00000090381',
        ndraws=100000,
        whichRank='rankeq',
        minScore=900,
        exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'STITCHslc6a17.pdf'),
        width=100,
        height=100)

ggBarcode(vdbr=vdbr,
          namesPath=here('annotateDrugDb', 'compounds.csv'),
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000090381',
          minScore=900,
          barwidth1=2,
          barwidth2=25,
          exportPath=here('predPharma', 'sorl1Enrich_f0exp1', 'bc_STITCHslc6a17.pdf'),
          width=200,
          height=100)


# human STITCH targets ----------------------------------------------------

# not ran yet...
tic()
hstitch <- drugEnrichment(vdbr=vdbr,
                          namesPath=here('annotateDrugDb', 'compounds.csv'),
                          annotationPath=here('annotateDrugDb', 'humanSTITCH.csv'),
                          annotation='humanSTITCH',
                          whichRank='rankeq',
                          minScore=900,
                          minNex=10,
                          ndraws=100000,
                          alphaThr=0.05,
                          maxPval=0.2,
                          statsExport=here('predPharma', 'sorl1Enrich_f0exp1', 'sorl1hSTITCH_f0exp1.csv'))
toc()


# zebrafish GO ------------------------------------------------------------

# not ran yet...
tic()
zgo <- drugEnrichment(vdbr=vdbr,
                      namesPath=here('annotateDrugDb', 'compounds.csv'),
                      annotationPath=here('annotateDrugDb', 'zebrafishGO.csv'),
                      annotation='zebrafishGO',
                      whichRank='rankeq',
                      minScore=900,
                      minNex=10,
                      ndraws=100000,
                      alphaThr=0.05,
                      maxPval=0.2,
                      statsExport=here('predPharma', 'sorl1Enrich_f0exp1', 'sorl1zGO_f0exp1.csv'))
toc()


# human GO ----------------------------------------------------------------

# not ran yet...
tic()
zgo <- drugEnrichment(vdbr=vdbr,
                      namesPath=here('annotateDrugDb', 'compounds.csv'),
                      annotationPath=here('annotateDrugDb', 'humanGO.csv'),
                      annotation='humanGO',
                      whichRank='rankeq',
                      minScore=900,
                      minNex=500,
                      ndraws=1000,
                      alphaThr=0.05,
                      maxPval=0.2,
                      statsExport=here('predPharma', 'sorl1Enrich_f0exp1', 'sorl1hGO_f0exp1.csv'))
toc()


### stable exp1 ###
lfp <- legacyFingerprint(matPath=here('220316_sorl1Stable/legacyMiddur/220316_14/220316_14.mat'),
                         conGrp='wt',
                         treGrp='hom',
                         days=c(2,3),
                         nights=c(2,3))
vdb <- rankDrugDb(legacyFgp=lfp,
                  dbPath=here('annotateDrugDb', 'drugDb.csv'), 
                  metric='cosine')
write.xlsx(vdb, here('predPharma', 'sorl1Enrich_staexp1', 'sorl1Ranked_staexp1.xlsx'))

tic()
zstitch <- drugEnrichment(lfp=lfp,
                          vdb=vdb, 
                          namesPath=here('annotateDrugDb', 'drugAnnotations_TTD.xlsx'),
                          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
                          annotation='zebrafishSTITCH',
                          minScore=700,
                          minNex=3,
                          ndraws=100000,
                          alphaThr=0.05,
                          maxPval=0.2)
toc()
write.xlsx(zstitch, here('predPharma', 'sorl1Enrich_f0exp1', 'sorl1zSTITCH_staexp1.xlsx'))


### stable exp2 ###
lfp <- legacyFingerprint(matPath=here('220316_sorl1Stable/legacyMiddur/220316_15/220316_15.mat'),
                         conGrp='wt',
                         treGrp='hom',
                         days=c(2,3),
                         nights=c(2,3))
vdb <- rankDrugDb(legacyFgp=lfp,
                  dbPath=here('annotateDrugDb', 'drugDb.csv'), 
                  metric='cosine')
write.xlsx(vdb, here('predPharma', 'sorl1Enrich_staexp2', 'sorl1Ranked_staexp2.xlsx'))

tic()
zstitch <- drugEnrichment(lfp=lfp,
                          vdb=vdb, 
                          namesPath=here('annotateDrugDb', 'drugAnnotations_TTD.xlsx'),
                          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
                          annotation='zebrafishSTITCH',
                          minScore=700,
                          minNex=3,
                          ndraws=100000,
                          alphaThr=0.05,
                          maxPval=0.2)
toc()
write.xlsx(zstitch, here('predPharma', 'sorl1Enrich_f0exp1', 'sorl1zSTITCH_staexp2.xlsx'))


# -------------------------------------------------------------------------

# with sorl1 f0 exp2 ------------------------------------------------------

dir.create(here('predPharma', 'sorl1Enrich_f0exp2'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_f0exp2', 'indications'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_f0exp2', 'targets'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_f0exp2', 'targetBioclass'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_f0exp2', 'keggPaths'), showWarnings=FALSE)

# calculate legacy fingerprint
lfp <- legacyFingerprint(matPath='~/Dropbox/ZFAD/220531_SORL1/legacyMiddur/220531_15/220531_15.mat',
                         conGrp='scr',
                         treGrp='sorl1',
                         days=c(2,3),
                         nights=c(2,3))

# rank drug db in comparison to fingerprint
vdb <- rankDrugDb(legacyFgp=lfp,
                  dbPath='~/Dropbox/phd/drugDatabase/drugDb.csv', 
                  metric='cosine')
write.xlsx(vdb, here('predPharma', 'sorl1Enrich_f0exp2', 'sorl1Ranked_f0exp2.xlsx'))

# indications
inr <- drugEnrichment(lfp=lfp,
                      vdb=vdb, 
                      annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                      annotation='indication',
                      minNex=3,
                      ndraws=10000,
                      alphaThr=0.05,
                      onlyPlotSign=FALSE,
                      plotDir=here('predPharma', 'sorl1Enrich_f0exp2', 'indications'))
write.xlsx(inr, here('predPharma', 'sorl1Enrich_f0exp2', 'sorl1Indications_f0exp2.xlsx'))

# targets
tar <- drugEnrichment(lfp=lfp,
                      vdb=vdb, 
                      annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                      annotation='target',
                      minNex=3,
                      ndraws=10000,
                      alphaThr=0.05,
                      onlyPlotSign=TRUE,
                      plotDir=here('predPharma', 'sorl1Enrich_f0exp2', 'targets'))
write.xlsx(tar, here('predPharma', 'sorl1Enrich_f0exp2', 'sorl1Targets_f0exp2.xlsx'))

# target bioclass
tarbio <- drugEnrichment(lfp=lfp,
                         vdb=vdb, 
                         annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                         annotation='targetBioclass',
                         minNex=3,
                         ndraws=10000,
                         alphaThr=0.05,
                         onlyPlotSign=TRUE,
                         plotDir=here('predPharma', 'sorl1Enrich_f0exp2', 'targetBioclass'))
write.xlsx(tarbio, here('predPharma', 'sorl1Enrich_f0exp2', 'sorl1TargetBioclass_f0exp2.xlsx'))

# KEGG pathways
kegg <- drugEnrichment(lfp=lfp,
                       vdb=vdb, 
                       annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                       annotation='keggPathways',
                       minNex=3,
                       ndraws=10000,
                       alphaThr=0.05,
                       onlyPlotSign=TRUE,
                       plotDir=here('predPharma', 'sorl1Enrich_f0exp2', 'keggPaths'))
write.xlsx(kegg, here('predPharma', 'sorl1Enrich_f0exp2', 'sorl1KEGG_f0exp2.xlsx'))



# with sorl1 stable exp1 --------------------------------------------------

dir.create(here('predPharma', 'sorl1Enrich_staexp1'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_staexp1', 'indications'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_staexp1', 'targets'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_staexp1', 'targetBioclass'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_staexp1', 'keggPaths'), showWarnings=FALSE)

# calculate legacy fingerprint
lfp <- legacyFingerprint(matPath='~/Dropbox/ZFAD/220316_sorl1Stable/legacyMiddur/220316_14/220316_14.mat',
                         conGrp='wt',
                         treGrp='hom',
                         days=c(2,3),
                         nights=c(2,3))

# rank drug db in comparison to fingerprint
vdb <- rankDrugDb(legacyFgp=lfp,
                  dbPath='~/Dropbox/phd/drugDatabase/drugDb.csv', 
                  metric='cosine')
write.xlsx(vdb, here('predPharma', 'sorl1Enrich_staexp1', 'sorl1Ranked_staexp1.xlsx'))

# indications
inr <- drugEnrichment(lfp=lfp,
                      vdb=vdb, 
                      annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                      annotation='indication',
                      minNex=3,
                      ndraws=10000,
                      alphaThr=0.05,
                      onlyPlotSign=FALSE,
                      plotDir=here('predPharma', 'sorl1Enrich_staexp1', 'indications'))
write.xlsx(inr, here('predPharma', 'sorl1Enrich_staexp1', 'sorl1Indications_staexp1.xlsx'))

# targets
tar <- drugEnrichment(lfp=lfp,
                      vdb=vdb, 
                      annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                      annotation='target',
                      minNex=3,
                      ndraws=10000,
                      alphaThr=0.05,
                      onlyPlotSign=TRUE,
                      plotDir=here('predPharma', 'sorl1Enrich_staexp1', 'targets'))
write.xlsx(tar, here('predPharma', 'sorl1Enrich_staexp1', 'sorl1Targets_staexp1.xlsx'))

# target bioclass
tarbio <- drugEnrichment(lfp=lfp,
                         vdb=vdb, 
                         annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                         annotation='targetBioclass',
                         minNex=3,
                         ndraws=10000,
                         alphaThr=0.05,
                         onlyPlotSign=TRUE,
                         plotDir=here('predPharma', 'sorl1Enrich_staexp1', 'targetBioclass'))
write.xlsx(tarbio, here('predPharma', 'sorl1Enrich_staexp1', 'sorl1TargetBioclass_staexp1.xlsx'))

# KEGG pathways
kegg <- drugEnrichment(lfp=lfp,
                       vdb=vdb, 
                       annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                       annotation='keggPathways',
                       minNex=3,
                       ndraws=10000,
                       alphaThr=0.05,
                       onlyPlotSign=TRUE,
                       plotDir=here('predPharma', 'sorl1Enrich_staexp1', 'keggPaths'))
write.xlsx(kegg, here('predPharma', 'sorl1Enrich_staexp1', 'sorl1KEGG_staexp1.xlsx'))


# with sorl1 stable exp2 --------------------------------------------------

dir.create(here('predPharma', 'sorl1Enrich_staexp2'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_staexp2', 'indications'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_staexp2', 'targets'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_staexp2', 'targetBioclass'), showWarnings=FALSE)
dir.create(here('predPharma', 'sorl1Enrich_staexp2', 'keggPaths'), showWarnings=FALSE)

# calculate legacy fingerprint
lfp <- legacyFingerprint(matPath=here('220316_sorl1Stable/legacyMiddur/220316_15/220316_15.mat'),
                         conGrp='wt',
                         treGrp='hom',
                         days=c(2,3),
                         nights=c(2,3))

# rank drug db in comparison to fingerprint
vdb <- rankDrugDb(legacyFgp=lfp,
                  dbPath=here('annotateDrugDb/drugDb.csv'), 
                  metric='cosine')
write.xlsx(vdb, here('predPharma', 'sorl1Enrich_staexp2', 'sorl1Ranked_staexp2.xlsx'))

# indications
inr <- drugEnrichment(lfp=lfp,
                      vdb=vdb, 
                      annotationPath=here('annotateDrugDb/drugAnnotations_TTD.xlsx'),
                      annotation='indication',
                      minNex=3,
                      ndraws=10000,
                      alphaThr=0.05)
write.xlsx(inr, here('predPharma', 'sorl1Enrich_staexp2', 'sorl1Indications_staexp2.xlsx'))


# targets
tar <- drugEnrichment(lfp=lfp,
                      vdb=vdb, 
                      annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                      annotation='target',
                      minNex=3,
                      ndraws=10000,
                      alphaThr=0.05,
                      onlyPlotSign=TRUE,
                      plotDir=here('predPharma', 'sorl1Enrich_staexp2', 'targets'))
write.xlsx(tar, here('predPharma', 'sorl1Enrich_staexp2', 'sorl1Targets_staexp2.xlsx'))

# target bioclass
tarbio <- drugEnrichment(lfp=lfp,
                         vdb=vdb, 
                         annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                         annotation='targetBioclass',
                         minNex=3,
                         ndraws=10000,
                         alphaThr=0.05,
                         onlyPlotSign=TRUE,
                         plotDir=here('predPharma', 'sorl1Enrich_staexp2', 'targetBioclass'))
write.xlsx(tarbio, here('predPharma', 'sorl1Enrich_staexp2', 'sorl1TargetBioclass_staexp2.xlsx'))

# KEGG pathways
kegg <- drugEnrichment(lfp=lfp,
                       vdb=vdb, 
                       annotationPath='~/Dropbox/phd/drugDatabase/drugAnnotations.xlsx',
                       annotation='keggPathways',
                       minNex=3,
                       ndraws=10000,
                       alphaThr=0.05,
                       onlyPlotSign=TRUE,
                       plotDir=here('predPharma', 'sorl1Enrich_staexp2', 'keggPaths'))
write.xlsx(kegg, here('predPharma', 'sorl1Enrich_staexp2', 'sorl1KEGG_staexp2.xlsx'))

