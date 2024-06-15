# 220524_miseq

Please follow in parallel 220524_MiSeq.main.R.  

Note; also copied data from separate September 2021 run, only csnk1db samples.  

See 220524_MiSeq_map1.xlsx & 220524_MiSeq_map2.xlsx for full source of samples.  

### Gather all fastq files

See 220524_MiSeq_main.R.

### Prepare fasta references

Amplicon sequences in zfad_seqs, column _amplicon_seq_.  

Fasta references are provided here in /fastarefs/ for reference but `prepareFastaRef.command` expects them in ~/Dropbox/phd/fastarefs.  

Each reference is 5'–3' amplicon, all lowercase except PAM sequence in uppercase.  

    prepareFastaRef.command apoea_1.fa
    prepareFastaRef.command apoea_2.fa

    prepareFastaRef.command apoeb_1.fa
    prepareFastaRef.command apoeb_3.fa

    prepareFastaRef.command appa_3.fa
    prepareFastaRef.command appa_4.fa

    prepareFastaRef.command appb_1.fa
    prepareFastaRef.command appb_2.fa

    prepareFastaRef.command sorl1_1.fa
    prepareFastaRef.command sorl1_2.fa
    prepareFastaRef.command sorl1_3.fa

    prepareFastaRef.command cd2ap_1.fa
    prepareFastaRef.command cd2ap_2.fa
    prepareFastaRef.command cd2ap_3.fa

    prepareFastaRef.command clu_1.fa
    prepareFastaRef.command clu_2.fa
    prepareFastaRef.command clu_3.fa

### Align reads to reference + Filter alignment file + Convert back to filtered fastq files

Filtering parameters are given within alignFilterBackSingle.command when it calls filterBam.command, as in: `filterBam.command -i $BAM -e 40 -f 100 -s 0.2 -p no -o $FILTER`

Filtering parameters were:

* minimum Phred score = 40 (i.e. delete any read with lower quality score)
* minimum reference span = 100 bp (i.e. delete any read that covers less of the reference)
* minimum proportion soft-clipped = 0.2 (i.e. delete any read that has a larger proportion soft-clipped)

Please go to scripts filterBam.command & readToThrow.R for complete explanations.  

Indeed, alignFilterBackSingle.command has dependencies I wrote, namely it calls filterBam.command, which calls readsToThrow.R.  

Copy-pasted fastq files we will use to new folder reads/.  

Commands:

    cd ~/Dropbox/phd/220524_miseq/data/reads/

    ### plate1

    alignFilterBackSingle.command B02_S14_L001_R1_001.fastq.gz apoea_1.fa
    alignFilterBackSingle.command B03_S15_L001_R1_001.fastq.gz apoea_1.fa
    alignFilterBackSingle.command B06_S18_L001_R1_001.fastq.gz apoea_1.fa
    alignFilterBackSingle.command B07_S19_L001_R1_001.fastq.gz apoea_1.fa
    alignFilterBackSingle.command B08_S20_L001_R1_001.fastq.gz apoea_1.fa
    alignFilterBackSingle.command B09_S21_L001_R1_001.fastq.gz apoea_1.fa
    alignFilterBackSingle.command B10_S22_L001_R1_001.fastq.gz apoea_1.fa
    alignFilterBackSingle.command B11_S23_L001_R1_001.fastq.gz apoea_1.fa

    alignFilterBackSingle.command C02_S26_L001_R1_001.fastq.gz apoea_2.fa
    alignFilterBackSingle.command C03_S27_L001_R1_001.fastq.gz apoea_2.fa
    alignFilterBackSingle.command C06_S30_L001_R1_001.fastq.gz apoea_2.fa
    alignFilterBackSingle.command C07_S31_L001_R1_001.fastq.gz apoea_2.fa
    alignFilterBackSingle.command C08_S32_L001_R1_001.fastq.gz apoea_2.fa
    alignFilterBackSingle.command C09_S33_L001_R1_001.fastq.gz apoea_2.fa
    alignFilterBackSingle.command C10_S34_L001_R1_001.fastq.gz apoea_2.fa
    alignFilterBackSingle.command C11_S35_L001_R1_001.fastq.gz apoea_2.fa

    alignFilterBackSingle.command D02_S38_L001_R1_001.fastq.gz apoeb_1.fa
    alignFilterBackSingle.command D03_S39_L001_R1_001.fastq.gz apoeb_1.fa
    alignFilterBackSingle.command D06_S42_L001_R1_001.fastq.gz apoeb_1.fa
    alignFilterBackSingle.command D07_S43_L001_R1_001.fastq.gz apoeb_1.fa
    alignFilterBackSingle.command D08_S44_L001_R1_001.fastq.gz apoeb_1.fa
    alignFilterBackSingle.command D09_S45_L001_R1_001.fastq.gz apoeb_1.fa
    alignFilterBackSingle.command D10_S46_L001_R1_001.fastq.gz apoeb_1.fa
    alignFilterBackSingle.command D11_S47_L001_R1_001.fastq.gz apoeb_1.fa

    alignFilterBackSingle.command E02_S50_L001_R1_001.fastq.gz apoeb_3.fa
    alignFilterBackSingle.command E03_S51_L001_R1_001.fastq.gz apoeb_3.fa
    alignFilterBackSingle.command E06_S54_L001_R1_001.fastq.gz apoeb_3.fa
    alignFilterBackSingle.command E07_S55_L001_R1_001.fastq.gz apoeb_3.fa
    alignFilterBackSingle.command E08_S56_L001_R1_001.fastq.gz apoeb_3.fa
    alignFilterBackSingle.command E09_S57_L001_R1_001.fastq.gz apoeb_3.fa
    alignFilterBackSingle.command E10_S58_L001_R1_001.fastq.gz apoeb_3.fa
    alignFilterBackSingle.command E11_S59_L001_R1_001.fastq.gz apoeb_3.fa

    alignFilterBackSingle.command F02_S62_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command F03_S63_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command F06_S66_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command F07_S67_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command F08_S68_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command F09_S69_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command F10_S70_L001_R1_001.fastq.gz appa_3.fa
    alignFilterBackSingle.command F11_S71_L001_R1_001.fastq.gz appa_3.fa

    alignFilterBackSingle.command G02_S74_L001_R1_001.fastq.gz appa_4.fa
    alignFilterBackSingle.command G03_S75_L001_R1_001.fastq.gz appa_4.fa
    alignFilterBackSingle.command G06_S78_L001_R1_001.fastq.gz appa_4.fa
    alignFilterBackSingle.command G07_S79_L001_R1_001.fastq.gz appa_4.fa
    alignFilterBackSingle.command G08_S80_L001_R1_001.fastq.gz appa_4.fa
    alignFilterBackSingle.command G09_S81_L001_R1_001.fastq.gz appa_4.fa
    alignFilterBackSingle.command G10_S82_L001_R1_001.fastq.gz appa_4.fa
    alignFilterBackSingle.command G11_S83_L001_R1_001.fastq.gz appa_4.fa

    ### plate2

    alignFilterBackSingle.command B02_S14_L001_R1_001.fastq.gz appb_1.fa
    alignFilterBackSingle.command B03_S15_L001_R1_001.fastq.gz appb_1.fa
    alignFilterBackSingle.command B06_S18_L001_R1_001.fastq.gz appb_1.fa
    alignFilterBackSingle.command B07_S19_L001_R1_001.fastq.gz appb_1.fa
    alignFilterBackSingle.command B08_S20_L001_R1_001.fastq.gz appb_1.fa
    alignFilterBackSingle.command B09_S21_L001_R1_001.fastq.gz appb_1.fa
    alignFilterBackSingle.command B10_S22_L001_R1_001.fastq.gz appb_1.fa
    alignFilterBackSingle.command B11_S23_L001_R1_001.fastq.gz appb_1.fa

    alignFilterBackSingle.command C02_S26_L001_R1_001.fastq.gz appb_2.fa
    alignFilterBackSingle.command C03_S27_L001_R1_001.fastq.gz appb_2.fa
    alignFilterBackSingle.command C06_S30_L001_R1_001.fastq.gz appb_2.fa
    alignFilterBackSingle.command C07_S31_L001_R1_001.fastq.gz appb_2.fa
    alignFilterBackSingle.command C08_S32_L001_R1_001.fastq.gz appb_2.fa
    alignFilterBackSingle.command C09_S33_L001_R1_001.fastq.gz appb_2.fa
    alignFilterBackSingle.command C10_S34_L001_R1_001.fastq.gz appb_2.fa
    alignFilterBackSingle.command C11_S35_L001_R1_001.fastq.gz appb_2.fa

    alignFilterBackSingle.command D02_S38_L001_R1_001.fastq.gz sorl1_1.fa
    alignFilterBackSingle.command D03_S39_L001_R1_001.fastq.gz sorl1_1.fa
    alignFilterBackSingle.command D06_S42_L001_R1_001.fastq.gz sorl1_1.fa
    alignFilterBackSingle.command D07_S43_L001_R1_001.fastq.gz sorl1_1.fa
    alignFilterBackSingle.command D08_S44_L001_R1_001.fastq.gz sorl1_1.fa
    alignFilterBackSingle.command D09_S45_L001_R1_001.fastq.gz sorl1_1.fa
    alignFilterBackSingle.command D10_S46_L001_R1_001.fastq.gz sorl1_1.fa
    alignFilterBackSingle.command D11_S47_L001_R1_001.fastq.gz sorl1_1.fa

    alignFilterBackSingle.command E02_S50_L001_R1_001.fastq.gz sorl1_2.fa
    alignFilterBackSingle.command E03_S51_L001_R1_001.fastq.gz sorl1_2.fa
    alignFilterBackSingle.command E06_S54_L001_R1_001.fastq.gz sorl1_2.fa
    alignFilterBackSingle.command E07_S55_L001_R1_001.fastq.gz sorl1_2.fa
    alignFilterBackSingle.command E08_S56_L001_R1_001.fastq.gz sorl1_2.fa
    alignFilterBackSingle.command E09_S57_L001_R1_001.fastq.gz sorl1_2.fa
    alignFilterBackSingle.command E10_S58_L001_R1_001.fastq.gz sorl1_2.fa
    alignFilterBackSingle.command E11_S59_L001_R1_001.fastq.gz sorl1_2.fa

    alignFilterBackSingle.command F02_S62_L001_R1_001.fastq.gz sorl1_3.fa
    alignFilterBackSingle.command F03_S63_L001_R1_001.fastq.gz sorl1_3.fa
    alignFilterBackSingle.command F06_S66_L001_R1_001.fastq.gz sorl1_3.fa
    alignFilterBackSingle.command F07_S67_L001_R1_001.fastq.gz sorl1_3.fa
    alignFilterBackSingle.command F08_S68_L001_R1_001.fastq.gz sorl1_3.fa
    alignFilterBackSingle.command F09_S69_L001_R1_001.fastq.gz sorl1_3.fa
    alignFilterBackSingle.command F10_S70_L001_R1_001.fastq.gz sorl1_3.fa
    alignFilterBackSingle.command F11_S71_L001_R1_001.fastq.gz sorl1_3.fa

    ### plate3

    alignFilterBackSingle.command B02_S14_L001_R1_001.fastq.gz cd2ap_1.fa
    alignFilterBackSingle.command B03_S15_L001_R1_001.fastq.gz cd2ap_1.fa
    alignFilterBackSingle.command B06_S18_L001_R1_001.fastq.gz cd2ap_1.fa
    alignFilterBackSingle.command B07_S19_L001_R1_001.fastq.gz cd2ap_1.fa
    alignFilterBackSingle.command B08_S20_L001_R1_001.fastq.gz cd2ap_1.fa
    alignFilterBackSingle.command B09_S21_L001_R1_001.fastq.gz cd2ap_1.fa
    alignFilterBackSingle.command B10_S22_L001_R1_001.fastq.gz cd2ap_1.fa
    alignFilterBackSingle.command B11_S23_L001_R1_001.fastq.gz cd2ap_1.fa

    alignFilterBackSingle.command C02_S26_L001_R1_001.fastq.gz cd2ap_2.fa
    alignFilterBackSingle.command C03_S27_L001_R1_001.fastq.gz cd2ap_2.fa
    alignFilterBackSingle.command C06_S30_L001_R1_001.fastq.gz cd2ap_2.fa
    alignFilterBackSingle.command C07_S31_L001_R1_001.fastq.gz cd2ap_2.fa
    alignFilterBackSingle.command C08_S32_L001_R1_001.fastq.gz cd2ap_2.fa
    alignFilterBackSingle.command C09_S33_L001_R1_001.fastq.gz cd2ap_2.fa
    alignFilterBackSingle.command C10_S34_L001_R1_001.fastq.gz cd2ap_2.fa
    alignFilterBackSingle.command C11_S35_L001_R1_001.fastq.gz cd2ap_2.fa

    alignFilterBackSingle.command D02_S38_L001_R1_001.fastq.gz cd2ap_3.fa
    alignFilterBackSingle.command D03_S39_L001_R1_001.fastq.gz cd2ap_3.fa
    alignFilterBackSingle.command D06_S42_L001_R1_001.fastq.gz cd2ap_3.fa
    alignFilterBackSingle.command D07_S43_L001_R1_001.fastq.gz cd2ap_3.fa
    alignFilterBackSingle.command D08_S44_L001_R1_001.fastq.gz cd2ap_3.fa
    alignFilterBackSingle.command D09_S45_L001_R1_001.fastq.gz cd2ap_3.fa
    alignFilterBackSingle.command D10_S46_L001_R1_001.fastq.gz cd2ap_3.fa
    alignFilterBackSingle.command D11_S47_L001_R1_001.fastq.gz cd2ap_3.fa

    alignFilterBackSingle.command E02_S50_L001_R1_001.fastq.gz clu_1.fa
    alignFilterBackSingle.command E03_S51_L001_R1_001.fastq.gz clu_1.fa
    alignFilterBackSingle.command E06_S54_L001_R1_001.fastq.gz clu_1.fa
    alignFilterBackSingle.command E07_S55_L001_R1_001.fastq.gz clu_1.fa
    alignFilterBackSingle.command E08_S56_L001_R1_001.fastq.gz clu_1.fa
    alignFilterBackSingle.command E09_S57_L001_R1_001.fastq.gz clu_1.fa
    alignFilterBackSingle.command E10_S58_L001_R1_001.fastq.gz clu_1.fa
    alignFilterBackSingle.command E11_S59_L001_R1_001.fastq.gz clu_1.fa

    alignFilterBackSingle.command F02_S62_L001_R1_001.fastq.gz clu_2.fa
    alignFilterBackSingle.command F03_S63_L001_R1_001.fastq.gz clu_2.fa
    alignFilterBackSingle.command F06_S66_L001_R1_001.fastq.gz clu_2.fa
    alignFilterBackSingle.command F07_S67_L001_R1_001.fastq.gz clu_2.fa
    alignFilterBackSingle.command F08_S68_L001_R1_001.fastq.gz clu_2.fa
    alignFilterBackSingle.command F09_S69_L001_R1_001.fastq.gz clu_2.fa
    alignFilterBackSingle.command F10_S70_L001_R1_001.fastq.gz clu_2.fa
    alignFilterBackSingle.command F11_S71_L001_R1_001.fastq.gz clu_2.fa

    alignFilterBackSingle.command G02_S74_L001_R1_001.fastq.gz clu_3.fa
    alignFilterBackSingle.command G03_S75_L001_R1_001.fastq.gz clu_3.fa
    alignFilterBackSingle.command G06_S78_L001_R1_001.fastq.gz clu_3.fa
    alignFilterBackSingle.command G07_S79_L001_R1_001.fastq.gz clu_3.fa
    alignFilterBackSingle.command G08_S80_L001_R1_001.fastq.gz clu_3.fa
    alignFilterBackSingle.command G09_S81_L001_R1_001.fastq.gz clu_3.fa
    alignFilterBackSingle.command G10_S82_L001_R1_001.fastq.gz clu_3.fa
    alignFilterBackSingle.command G11_S83_L001_R1_001.fastq.gz clu_3.fa




Note; E06.bam vs apoeb_3 has a bunch of mismatches I think created by Cas9, might not be counted by ampliCan.
appa_4: contamination in controls. They have mutations.
Also weird repair. Mutations are two sites. Make sure ampliCan gets both.

Can copy-paste all above commands in script `batchCommands.command` and run it from Terminal.   

### ampliCan

See 220525_MiSeq.main.R.  

Note, turned off normalisation to controls by not labelling control samples in config file.  Rationale: it can decrease mutation rate because it finds mutations in control samples which is clearly contamination, while it is mostly useless: we are plotting only indels and it is extremely rare that there is an indel naturally present in the buffer around the Cas9 double-strand break site. Even if that was the case would detect it if it was the case and then we could take action.  

#### Inspecting samples with problematic coverage

##### plate 1

* apoea.1_scr1 = B03: no coverage before filtering  


* apoeb.1_ko5 = D10: two bad reads before filtering


* apoeb.3_ko3 = E08: low coverage before filtering
* apoeb.3_ko5 = E10: no reads before filtering

##### plate 2

* sorl1.3_ko4 = F9: no coverage before filtering
* sorl1.3_ko5 = F10: no coverage before filtering

##### plate 3

* cd2ap.1_scr1 = B02: no reads before filtering, band on gel seemed ok  
* cd2ap.1_scr2 = B03: lots of clipped reads, 28 pairs left after filtering
* cd2ap.1_ko2 = B07: low coverage already before filtering
* cd2ap.1_ko3 = B08: no reads already before filtering
* cd2ap.1_ko5 = B10: no reads already before filtering
* cd2ap.1_ko6 = B11: only 1 pair of reads before filtering  
Are ok: cd2ap.1_scr2 and cd2ap.1_ko1 and cd2ap.1_ko4 samples. Can send new samples but should be enough.  


* clu.1_scr1 = E02
* clu.1_scr2 = E03
* clu.1_ko1 = E06
* clu.1_ko2 = E07
* clu.1_ko3 = E08
* clu.1_ko4 = E09
* clu.1_ko5 = E10
* clu.1_ko6 = E11  
All low or no coverage already before filtering.  
Reference sequence looks correct.  
There was no band on gel! Will have to design new primers and send again.  


* clu.3_scr1 = G2: already low coverage before filtering
* clu.3_scr2 = G3: no coverage before filtering
* clu.3_ko3 = G8: no coverage before filtering
* clu.3_ko4 = G9: no coverage before filtering
* clu.3_ko5 = G10: **
* clu.3_ko6 = G11: **

Note, gel did not look good for clu.3_scr.  
Ideally would need to send new samples.   

\** Lots of trimmed reads before filtering. Stops matching exactly where expected double-strand break – is it a big deletion that removed Forward primer? Soft-clipped portion is  
ACAGCGTTCTGTTTGTAGCGATCCAAAGCT,
which is Forward primer (ACAGCGTTCTGTTTGTAGCGA) + 9 bp of start of the amplicon. Read https://www.biostars.org/p/9529235/   
