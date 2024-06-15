# October 2021 MiSeq

Please follow in parallel oct2021_MiSeq.main.R.  

Note; also copied data from separate September 2021 run, only csnk1db samples.  

All csnk1db samples are 210412_BOX11.

See octMiSeq_map1.xlsx & octMiSeq_map2.xlsx for full source of sample (can match with behaviour data).  


## About September 2021 MiSeq

Sent exactly same csnk1db amplicons as used in headloop on 02/09/2021; except
* headloop PCR, as headloop primers did not have MiSeq tags  
* all larva KO_4, as it did not produce any band  

### Gather all fastq files

See oct2021_MiSeq.main.R.

### Prepare fasta references

Fasta references for _csnk1db_ targets were prepared during F0 knockout project.

Fasta references are provided here in oct2021_fastarefs/ for reference, but code expects them in ~/Dropbox/phd/fastarefs.  

Each reference is 5'–3' amplicon, all lowercase except PAM sequence in uppercase.  

    prepareFastaRef.command appa_1.fa
    prepareFastaRef.command appa_2.fa
    prepareFastaRef.command appa_3.fa

    prepareFastaRef.command psen1_1.fa
    prepareFastaRef.command psen1_2.fa
    prepareFastaRef.command psen1_3.fa

    prepareFastaRef.command psen2_1.fa
    prepareFastaRef.command psen2_2.fa
    prepareFastaRef.command psen2_3.fa

    prepareFastaRef.command csnk1db_a.fa
    prepareFastaRef.command csnk1db_b.fa
    prepareFastaRef.command csnk1db_d.fa

### Align reads to reference + Filter alignment file + Convert back to filtered fastq files

Filtering parameters are given within alignFilterBackSingle.command when it calls filterBam.command, as in: `filterBam.command -i $BAM -e 40 -f 100 -s 0.2 -p no -o $FILTER`

Filtering parameters were:

* minimum Phred score = 40 (i.e. delete any read with lower quality score)
* minimum reference span = 100 bp (i.e. delete any read that covers less of the reference)
* minimum proportion soft-clipped = 0.2 (i.e. delete any read that has a larger proportion soft-clipped)

Please go to scripts filterBam.command & readToThrow.R for complete explanations.  

Indeed, alignFilterBackSingle.command has dependencies I wrote, namely it calls filterBam.command, which calls readsToThrow.R.  

#### September2021

Copy-pasted fastq files we will use to new folder reads/.  

Commands:

    cd ~/Dropbox/phd/october2021_MiSeq/sept2021_csnk1db/reads/

    alignFilterBackSingle.command A02_S2_L001_R1_001.fastq.gz ref_csnk1db_aa.fa
    alignFilterBackSingle.command A03_S3_L001_R1_001.fastq.gz ref_csnk1db_aa.fa
    alignFilterBackSingle.command A04_S4_L001_R1_001.fastq.gz ref_csnk1db_aa.fa
    alignFilterBackSingle.command A05_S5_L001_R1_001.fastq.gz ref_csnk1db_aa.fa
    alignFilterBackSingle.command A06_S6_L001_R1_001.fastq.gz ref_csnk1db_aa.fa

    alignFilterBackSingle.command B02_S14_L001_R1_001.fastq.gz ref_csnk1db_ab.fa
    alignFilterBackSingle.command B03_S15_L001_R1_001.fastq.gz ref_csnk1db_ab.fa
    alignFilterBackSingle.command B04_S16_L001_R1_001.fastq.gz ref_csnk1db_ab.fa
    alignFilterBackSingle.command B05_S17_L001_R1_001.fastq.gz ref_csnk1db_ab.fa
    alignFilterBackSingle.command B06_S18_L001_R1_001.fastq.gz ref_csnk1db_ab.fa

    alignFilterBackSingle.command C02_S26_L001_R1_001.fastq.gz ref_csnk1db_ad.fa
    alignFilterBackSingle.command C03_S27_L001_R1_001.fastq.gz ref_csnk1db_ad.fa
    alignFilterBackSingle.command C04_S28_L001_R1_001.fastq.gz ref_csnk1db_ad.fa
    alignFilterBackSingle.command C05_S29_L001_R1_001.fastq.gz ref_csnk1db_ad.fa
    alignFilterBackSingle.command C06_S30_L001_R1_001.fastq.gz ref_csnk1db_ad.fa

Can copy-paste all above commands in script `batchCommands.command` and run it from Terminal.   

#### October2021

Copy-pasted fastq files we will use to new folders plate1_data/ & plate2_data/.  

Note; some wells containing different amplicons were pooled for sequencing, so important to separate correctly in two folders plate1/plate2 below so files do not overwrite each other.

    ###

    cd ~/Dropbox/phd/october2021_MiSeq/plate1_data

    # psen1

    alignFilterBackSingle.command A02_S2_L001_R1_001.fastq.gz psen1_1.fa
    alignFilterBackSingle.command A03_S3_L001_R1_001.fastq.gz psen1_1.fa
    alignFilterBackSingle.command A04_S4_L001_R1_001.fastq.gz psen1_1.fa
    alignFilterBackSingle.command A06_S6_L001_R1_001.fastq.gz psen1_1.fa
    alignFilterBackSingle.command A07_S7_L001_R1_001.fastq.gz psen1_1.fa
    alignFilterBackSingle.command A08_S8_L001_R1_001.fastq.gz psen1_1.fa
    alignFilterBackSingle.command A09_S9_L001_R1_001.fastq.gz psen1_1.fa
    alignFilterBackSingle.command A10_S10_L001_R1_001.fastq.gz psen1_1.fa

    alignFilterBackSingle.command B02_S14_L001_R1_001.fastq.gz psen1_2.fa
    alignFilterBackSingle.command B03_S15_L001_R1_001.fastq.gz psen1_2.fa
    alignFilterBackSingle.command B04_S16_L001_R1_001.fastq.gz psen1_2.fa
    alignFilterBackSingle.command B06_S18_L001_R1_001.fastq.gz psen1_2.fa
    alignFilterBackSingle.command B07_S19_L001_R1_001.fastq.gz psen1_2.fa
    alignFilterBackSingle.command B08_S20_L001_R1_001.fastq.gz psen1_2.fa
    alignFilterBackSingle.command B09_S21_L001_R1_001.fastq.gz psen1_2.fa
    alignFilterBackSingle.command B10_S22_L001_R1_001.fastq.gz psen1_2.fa

    alignFilterBackSingle.command C02_S26_L001_R1_001.fastq.gz psen1_3.fa
    alignFilterBackSingle.command C03_S27_L001_R1_001.fastq.gz psen1_3.fa
    alignFilterBackSingle.command C04_S28_L001_R1_001.fastq.gz psen1_3.fa
    alignFilterBackSingle.command C06_S30_L001_R1_001.fastq.gz psen1_3.fa
    alignFilterBackSingle.command C07_S31_L001_R1_001.fastq.gz psen1_3.fa
    alignFilterBackSingle.command C08_S32_L001_R1_001.fastq.gz psen1_3.fa
    alignFilterBackSingle.command C09_S33_L001_R1_001.fastq.gz psen1_3.fa
    alignFilterBackSingle.command C10_S34_L001_R1_001.fastq.gz psen1_3.fa

    ###

    # psen2

    cd ~/Dropbox/phd/october2021_MiSeq/plate1_data

    alignFilterBackSingle.command D02_S38_L001_R1_001.fastq.gz psen2_1.fa
    alignFilterBackSingle.command D03_S39_L001_R1_001.fastq.gz psen2_1.fa
    alignFilterBackSingle.command D04_S40_L001_R1_001.fastq.gz psen2_1.fa
    alignFilterBackSingle.command D06_S42_L001_R1_001.fastq.gz psen2_1.fa
    alignFilterBackSingle.command D07_S43_L001_R1_001.fastq.gz psen2_1.fa
    alignFilterBackSingle.command D08_S44_L001_R1_001.fastq.gz psen2_1.fa
    alignFilterBackSingle.command D09_S45_L001_R1_001.fastq.gz psen2_1.fa
    alignFilterBackSingle.command D10_S46_L001_R1_001.fastq.gz psen2_1.fa

    alignFilterBackSingle.command E02_S50_L001_R1_001.fastq.gz psen2_2.fa
    alignFilterBackSingle.command E03_S51_L001_R1_001.fastq.gz psen2_2.fa
    alignFilterBackSingle.command E04_S52_L001_R1_001.fastq.gz psen2_2.fa
    alignFilterBackSingle.command E06_S54_L001_R1_001.fastq.gz psen2_2.fa
    alignFilterBackSingle.command E07_S55_L001_R1_001.fastq.gz psen2_2.fa
    alignFilterBackSingle.command E08_S56_L001_R1_001.fastq.gz psen2_2.fa
    alignFilterBackSingle.command E09_S57_L001_R1_001.fastq.gz psen2_2.fa
    alignFilterBackSingle.command E10_S58_L001_R1_001.fastq.gz psen2_2.fa

    alignFilterBackSingle.command F02_S62_L001_R1_001.fastq.gz psen2_3.fa
    alignFilterBackSingle.command F03_S63_L001_R1_001.fastq.gz psen2_3.fa
    alignFilterBackSingle.command F04_S64_L001_R1_001.fastq.gz psen2_3.fa
    alignFilterBackSingle.command F06_S66_L001_R1_001.fastq.gz psen2_3.fa
    alignFilterBackSingle.command F07_S67_L001_R1_001.fastq.gz psen2_3.fa
    alignFilterBackSingle.command F08_S68_L001_R1_001.fastq.gz psen2_3.fa
    alignFilterBackSingle.command F09_S69_L001_R1_001.fastq.gz psen2_3.fa
    alignFilterBackSingle.command F10_S70_L001_R1_001.fastq.gz psen2_3.fa

    ###

    # csnk1db

    cd ~/Dropbox/phd/october2021_MiSeq/plate1_data

    alignFilterBackSingle.command G06_S78_L001_R1_001.fastq.gz csnk1db_b.fa
    alignFilterBackSingle.command G07_S79_L001_R1_001.fastq.gz csnk1db_b.fa

    alignFilterBackSingle.command H06_S90_L001_R1_001.fastq.gz csnk1db_d.fa
    alignFilterBackSingle.command H07_S91_L001_R1_001.fastq.gz csnk1db_d.fa

    ###

    # appa

    cd ~/Dropbox/phd/october2021_MiSeq/plate2_data

    alignFilterBackSingle.command C02_S26_L001_R1_001.fastq.gz appa_1.fa
    alignFilterBackSingle.command C03_S27_L001_R1_001.fastq.gz appa_1.fa
    alignFilterBackSingle.command C04_S28_L001_R1_001.fastq.gz appa_1.fa
    alignFilterBackSingle.command C06_S30_L001_R1_001.fastq.gz appa_1.fa
    alignFilterBackSingle.command C07_S31_L001_R1_001.fastq.gz appa_1.fa
    alignFilterBackSingle.command C08_S32_L001_R1_001.fastq.gz appa_1.fa
    alignFilterBackSingle.command C09_S33_L001_R1_001.fastq.gz appa_1.fa
    alignFilterBackSingle.command C10_S34_L001_R1_001.fastq.gz appa_1.fa

    alignFilterBackSingle.command D02_S38_L001_R1_001.fastq.gz appa_2.fa
    alignFilterBackSingle.command D03_S39_L001_R1_001.fastq.gz appa_2.fa
    alignFilterBackSingle.command D04_S40_L001_R1_001.fastq.gz appa_2.fa
    alignFilterBackSingle.command D06_S42_L001_R1_001.fastq.gz appa_2.fa
    alignFilterBackSingle.command D07_S43_L001_R1_001.fastq.gz appa_2.fa
    alignFilterBackSingle.command D08_S44_L001_R1_001.fastq.gz appa_2.fa
    alignFilterBackSingle.command D09_S45_L001_R1_001.fastq.gz appa_2.fa
    alignFilterBackSingle.command D10_S46_L001_R1_001.fastq.gz appa_2.fa

    alignFilterBackSingle.command E02_S50_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command E03_S51_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command E04_S52_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command E06_S54_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command E09_S57_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command E10_S58_L001_R1_001.fastq.gz appa_3.fa

### ampliCan

See oct2021_MiSeq.main.R.  

Checking samples with low coverage:

Plate 1  

* psen2.1_scr_3: many reads excluded because of soft-clipping
* psen1.2_ko_5: low coverage already before filtering
* csnk1db_D_ko_4: empty already before filtering
* csnk1db_B_ko_3: empty already before filtering

Plate 2  

* appa.3_ko_8: empty already before filtering; looked OK on gel though  
* appa.3_ko_7: empty already before filtering
