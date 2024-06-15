# 220915_miseq

Please follow in parallel 220915_MiSeq.main.R.  

See 220915_MiSeq_map1.xlsx for full source of samples.  

### Gather all fastq files

See 220915_MiSeq_main.R.

### Prepare fasta references

Amplicon sequences in zfad_seqs, column _amplicon_seq_.  

Fasta references are provided here in /fastarefs/ for reference but `prepareFastaRef.command` expects them in ~/Dropbox/phd/fastarefs.  

Each reference is 5'–3' amplicon, all lowercase except PAM sequence in uppercase.  

    prepareFastaRef.command ref_apoeb_ae.fa

    prepareFastaRef.command clu_1v2.fa
    prepareFastaRef.command clu_3v2.fa

    prepareFastaRef.command cd2ap_1v2.fa

### Align reads to reference + Filter alignment file + Convert back to filtered fastq files

Filtering parameters are given within alignFilterBackSingle.command when it calls filterBam.command, currently command is: `filterBam.command -i $BAM -d 20 -p no -o $FILTER`

Simplified filtering so that only criterion is that read has to cover double-strand break position ± 20 bp.

Please go to scripts filterBam.command & readToThrow.R for complete explanations.  

Indeed, alignFilterBackSingle.command has dependencies I wrote, namely it calls filterBam.command, which calls readsToThrow.R.  

Copy-pasted fastq files we need to new folder reads/.  

Commands:

    cd ~/Dropbox/phd/220915_miseq/data/reads/

    alignFilterBackSingle.command B02_S14_L001_R1_001.fastq.gz clu_1v2.fa
    alignFilterBackSingle.command B03_S15_L001_R1_001.fastq.gz clu_1v2.fa
    alignFilterBackSingle.command B04_S16_L001_R1_001.fastq.gz clu_1v2.fa
    alignFilterBackSingle.command B06_S18_L001_R1_001.fastq.gz clu_1v2.fa
    alignFilterBackSingle.command B07_S19_L001_R1_001.fastq.gz clu_1v2.fa
    alignFilterBackSingle.command B08_S20_L001_R1_001.fastq.gz clu_1v2.fa
    alignFilterBackSingle.command B09_S21_L001_R1_001.fastq.gz clu_1v2.fa
    alignFilterBackSingle.command B10_S22_L001_R1_001.fastq.gz clu_1v2.fa
    alignFilterBackSingle.command B11_S23_L001_R1_001.fastq.gz clu_1v2.fa

    alignFilterBackSingle.command C02_S26_L001_R1_001.fastq.gz clu_3v2.fa
    alignFilterBackSingle.command C03_S27_L001_R1_001.fastq.gz clu_3v2.fa
    alignFilterBackSingle.command C04_S28_L001_R1_001.fastq.gz clu_3v2.fa
    alignFilterBackSingle.command C06_S30_L001_R1_001.fastq.gz clu_3v2.fa
    alignFilterBackSingle.command C07_S31_L001_R1_001.fastq.gz clu_3v2.fa
    alignFilterBackSingle.command C08_S32_L001_R1_001.fastq.gz clu_3v2.fa
    alignFilterBackSingle.command C09_S33_L001_R1_001.fastq.gz clu_3v2.fa
    alignFilterBackSingle.command C10_S34_L001_R1_001.fastq.gz clu_3v2.fa
    alignFilterBackSingle.command C11_S35_L001_R1_001.fastq.gz clu_3v2.fa

    alignFilterBackSingle.command E02_S50_L001_R1_001.fastq.gz cd2ap_1v2.fa
    alignFilterBackSingle.command E03_S51_L001_R1_001.fastq.gz cd2ap_1v2.fa
    alignFilterBackSingle.command E04_S52_L001_R1_001.fastq.gz cd2ap_1v2.fa
    alignFilterBackSingle.command E06_S54_L001_R1_001.fastq.gz cd2ap_1v2.fa
    alignFilterBackSingle.command E07_S55_L001_R1_001.fastq.gz cd2ap_1v2.fa
    alignFilterBackSingle.command E08_S56_L001_R1_001.fastq.gz cd2ap_1v2.fa
    alignFilterBackSingle.command E09_S57_L001_R1_001.fastq.gz cd2ap_1v2.fa
    alignFilterBackSingle.command E10_S58_L001_R1_001.fastq.gz cd2ap_1v2.fa
    alignFilterBackSingle.command E11_S59_L001_R1_001.fastq.gz cd2ap_1v2.fa

    alignFilterBackSingle.command G02_S74_L001_R1_001.fastq.gz ref_apoeb_ae.fa
    alignFilterBackSingle.command G03_S75_L001_R1_001.fastq.gz ref_apoeb_ae.fa
    alignFilterBackSingle.command G04_S76_L001_R1_001.fastq.gz ref_apoeb_ae.fa
    alignFilterBackSingle.command G05_S77_L001_R1_001.fastq.gz ref_apoeb_ae.fa
    alignFilterBackSingle.command G06_S78_L001_R1_001.fastq.gz ref_apoeb_ae.fa
    alignFilterBackSingle.command G07_S79_L001_R1_001.fastq.gz ref_apoeb_ae.fa
    alignFilterBackSingle.command G08_S80_L001_R1_001.fastq.gz ref_apoeb_ae.fa
    alignFilterBackSingle.command G09_S81_L001_R1_001.fastq.gz ref_apoeb_ae.fa
    alignFilterBackSingle.command G10_S82_L001_R1_001.fastq.gz ref_apoeb_ae.fa
    alignFilterBackSingle.command G11_S83_L001_R1_001.fastq.gz ref_apoeb_ae.fa

Can copy-paste all above commands in script `batchCommands.command` and run it from Terminal.   

### ampliCan

See 220915_MiSeq.main.R.  

Note, turned off normalisation to controls by not labelling control samples in config file.  Rationale: it can decrease mutation rate because it finds mutations in control samples which is clearly contamination, while it is mostly useless: we are plotting only indels and it is extremely rare that there is an indel naturally present in the buffer around the Cas9 double-strand break site. Even if that was the case would detect it if it was the case and then we could take action.  

#### about apoeb samples

Question was whether we still had a SNP in background at chr16:23,961,572.  
Answer: we do not. There may still be carrier fish in the tank but unlikely.  
Can assume clean 10-bp deletion.  
KASP ordered with Eirinn's tool https://kasp.eirinn.org/ was designed with this assumption.  
