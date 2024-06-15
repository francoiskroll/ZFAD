#####################################################
# ~ ZFAD: generate transcript map for every gene tested ~
#
# using function transcriptMapper (in utilities/)
# Francois Kroll 2023
# francois@kroll.be
#####################################################



# source ------------------------------------------------------------------

library(here)
source(here('utilities', 'transcriptMapper_v3.R'))

dir.create(here('transcriptMaps', 'stables'), showWarnings=FALSE)



# sorl1 -------------------------------------------------------------------

transcriptMapper(transcript_id='ENSDART00000156995',
                 positions_tomark=c(21669448, 21594155, 21561866),
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=0.5,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps'),
                 width=70,
                 height=10)


# psen1 -------------------------------------------------------------------

# using psen1-201, gRNAs were designed with this one as reference

transcriptMapper(transcript_id='ENSDART00000149864',
                 positions_tomark=c(51218545, 51216830, 51210589),
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=1,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps'),
                 width=70,
                 height=10)

# psen2 -------------------------------------------------------------------

# using psen2-201, gRNAs were designed with this one as reference
transcriptMapper(transcript_id='ENSDART00000006381',
                 positions_tomark=c(51058244, 51056797, 51055741),
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=1.2,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps'),
                 width=70,
                 height=10)



# clu ---------------------------------------------------------------------

# using clu-201, gRNAs were designed with this one as reference
transcriptMapper(transcript_id='ENSDART00000127173',
                 positions_tomark=c(39271754, 39269439, 39263244),
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=1.2,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps'),
                 width=70,
                 height=10)



# clu stable --------------------------------------------------------------

# using clu-201, gRNAs were designed with this one as reference
transcriptMapper(transcript_id='ENSDART00000127173',
                 positions_tomark=c(39271603, 39271653),
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=1.2,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps', 'stables'),
                 width=70,
                 height=10)



# cd2ap -------------------------------------------------------------------

# using cd2ap-201, gRNAs were designed with this one as reference
transcriptMapper(transcript_id='ENSDART00000102611',
                 positions_tomark=c(35929963, 35936014, 35963036),
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=0,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps'),
                 width=70,
                 height=10)



# cd2ap stable ------------------------------------------------------------

# using cd2ap-201, gRNAs were designed with this one as reference
transcriptMapper(transcript_id='ENSDART00000102611',
                 positions_tomark=35900669,
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=0,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps', 'stables'),
                 width=70,
                 height=10)



# apoea -------------------------------------------------------------------

# using apoea-202, longest one
# for double F0 KO apoea/apoeb
# I used apoeb gRNA 1 & 3
transcriptMapper(transcript_id='ENSDART00000172219',
                 positions_tomark=c(10855890, 10856075),
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=1.2,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps'),
                 width=70,
                 height=10)



# apoeb -------------------------------------------------------------------

# using apoeb-201, longest one
# for double F0 KO apoea/apoeb
# I used apoeb gRNA 1 & 3
transcriptMapper(transcript_id='ENSDART00000058965',
                 positions_tomark=c(23962034, 23962861),
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=1.2,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps'),
                 width=70,
                 height=10)



# apoea stable ------------------------------------------------------------

# using apoea-202, longest one
transcriptMapper(transcript_id='ENSDART00000172219',
                 positions_tomark=10856109,
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=1.2,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps', 'stables'),
                 width=70,
                 height=10)


# apoeb stable ------------------------------------------------------------

# using apoeb-201, longest one
transcriptMapper(transcript_id='ENSDART00000058965',
                 positions_tomark=23961573,
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=1.2,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps', 'stables'),
                 width=70,
                 height=10)



# appa --------------------------------------------------------------------

# using appa-202
# in appa/b F0 experiments, used gRNA 3 & 4
transcriptMapper(transcript_id='ENSDART00000166786',
                 positions_tomark=c(625746, 633369),
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=0,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps'),
                 width=70,
                 height=10)



# appb --------------------------------------------------------------------

# using appb-203
# in appa/b F0 experiments, used gRNA 1 & 2
transcriptMapper(transcript_id='ENSDART00000077908',
                 positions_tomark=c(35059247, 35032746),
                 utr_height=5,
                 exon_height=8,
                 intron_height=1.5,
                 scale_exons=0.1,
                 scale_introns=0.01,
                 exon_colour='#697a87',
                 intron_colour='#aeb3b4',
                 roundCorners=TRUE,
                 corner_radius=0,
                 arrowsize=20,
                 arrowpos='below',
                 exportFolder=here('transcriptMaps'),
                 width=70,
                 height=10)
