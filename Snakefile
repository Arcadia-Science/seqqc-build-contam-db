import pandas as pd

genome_acc_df = pd.read_csv("inputs/genome_contaminants.csv") 
GENOME_ACC = genome_acc_df['ident'].unique().tolist()
sra_acc_df = pd.read_csv("inputs/seq_contaminants.csv")
SRA_ACC = sra_acc_df['ident'].unique().tolist()
KSIZE = [21]

rule all:
    input: expand("outputs/sourmash_contam_db/contamdb.dna.k{ksize}.zip", ksize = KSIZE)

#############################################
## SUBSET GTDB REPS DATABASE
## covers:
## - kit contamination
## - human microbiome contamination
#############################################

rule download_gtdb_rs207_metadata:
    output: "inputs/gtdb/bac120_metadata_r207.tar.gz"
    shell:'''
    curl -JLo {output} https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_metadata_r207.tar.gz
    '''

rule decompress_gtdb_rs207_metadata:
    input: "inputs/gtdb/bac120_metadata_r207.tar.gz"
    output: "inputs/gtdb/bac120_metadata_r207.tsv"
    params: outdir = "inputs/gtdb/"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule download_gtdb_rs207_taxonomy:
    output: "inputs/sourmash_databases/gtdb-rs207.taxonomy.csv.gz"
    shell:'''
    curl -JLo {output} https://osf.io/v3zmg/download
    '''

rule download_sourmash_db:
    input: "inputs/sourmash_databases/sourmash-database-info.csv"
    output: "inputs/sourmash_databases/gtdb-rs207.genomic-reps.dna.k{ksize}.zip"
    run:
        sourmash_database_info = pd.read_csv(str(input[0]))
        ksize = int(wildcards.ksize)
        db_df = sourmash_database_info.loc[(sourmash_database_info['ksize'] == ksize)]
        osf_hash = db_df['osf_hash'].values[0] 
        shell("curl -JLo {output} https://osf.io/{osf_hash}/download")


rule sourmash_sig_describe_db:
    '''
    produce a csv file that describes each of the signatures in the database.
    this CSV file will be post-processed to use as a picklist to extract signatures of interest form the sourmash db
    '''
    input: "inputs/sourmash_databases/gtdb-rs207.genomic-reps.dna.k{ksize}.zip"
    output: "outputs/sourmash_sig_describe/gtdb-rs207.genomic-reps.dna.k{ksize}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig describe -k {wildcards.ksize} --csv {output} {input}
    '''

rule subset_sourmash_db:
    input: 
        sig_describe = "outputs/sourmash_sig_describe/gtdb-rs207.genomic-reps.dna.k{ksize}.csv",
        contams = "inputs/doi10.1016j.tim.2018.11.003-table1.csv",
        gtdb_metadata = "inputs/gtdb/bac120_metadata_r207.tsv"
    output: picklist="outputs/sourmash_sig_describe_subset/picklist_k{ksize}.csv"
    conda: "envs/tidyverse.yml"
    script: "scripts/snakemake_subset_sourmash_db.R"

rule sourmash_extract_sigs:
    input: 
        db = "inputs/sourmash_databases/gtdb-rs207.genomic-reps.dna.k{ksize}.zip",
        picklist = "outputs/sourmash_sig_describe_subset/picklist_k{ksize}.csv"
    output: "outputs/sourmash_contam_db/gtdb-rs207.genomic-reps.dna.k{ksize}.contamsubset.zip"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig extract -o {output} --picklist {input.picklist}:md5:md5 {input.db}
    '''

##################################################################
## DOWNLOAD & SKETCH CONTAMINANT GENOMES
## covers:
## - commonly sequenced (model) species/index and barcode hopping
## - human contamination
## - some Arcadia-specific organisms
##################################################################

rule download_genomes:
    '''
    This rule downloads genomes that we expect might be common contaminants.
    Genomes are downloaded as zip files that contain the genome sequence, RNA seq (if available), and CDS from genomic (if available)
    '''
    output: "inputs/genomes/{genome_acc}.zip"
    shell:'''
    curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/{wildcards.genome_acc}/download?include_annotation_type=RNA_FASTA,CDS_FASTA&filename={wildcards.genome_acc}.zip" -H "Accept: application/zip"
    mv {wildcards.genome_acc}.zip inputs/genomes/{wildcards.genome_acc}.zip
    '''

rule unzip_genomes:
    '''
    Extract the sequence that will be sketched
    '''
    input: "inputs/genomes/{genome_acc}.zip"
    output: "inputs/genomes/{genome_acc}_cds_from_genomic.fna"
    conda: "envs/unzip.yml"
    shell:'''
    unzip -p {input} ncbi_dataset/data/{wildcards.genome_acc}/cds_from_genomic.fna > {output}
    '''

rule sourmash_sketch_genomes:
    '''
    Sketch the sequence.
    We settled on using the CDS from genomic sequence.
    Because eukaryotic genomes have a lot of repeat sequences, these might cause false positives when we search against them.
    Using RNA or CDS from genomic limits these false positives.
    All genomes were were interested in had CDS from genomic (but not RNA), so we settled on using this sequence.
    '''
    input: "inputs/genomes/{genome_acc}_cds_from_genomic.fna"
    output: "outputs/sourmash_sketch/{genome_acc}_cds_from_genomic.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund -o {output} --name {wildcards.genome_acc} {input}
    '''

##################################################################
## DOWNLOAD & SKETCH SEQUENCE READ ARCHIVE ACCESSIONS
## covers:
## - some Arcadia-specific organisms
##################################################################
    
rule sourmash_sketch_fastqs:
    '''
    some organisms that we expect to be contaminants don't have genomes on GenBank, but they do have raw seq data from genomes or transcriptomes.
    this rule downloads that data and generates a sketch from it.
    '''
    output: "outputs/sourmash_sketch/{sra_acc}_sra.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    fastq-dump --disable-multithreading --fasta 0 --skip-technical --readids --read-filter pass --dumpbase --split-spot --clip -Z {wildcards.sra_acc} | 
        sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {wildcards.sra_acc} -o {output} -
    '''

rule sourmash_filter:
    '''
    filter out hashes that only have abundance 1. 
    because these signatures are from raw fastq files, anything with abundance 1 is probably a sequencing error.
    removing these will decrease the database size and increase search times
    '''
    input: "outputs/sourmash_sketch/{sra_acc}_sra.sig"
    output: "outputs/sourmash_sketch_filtered/{sra_acc}_sra.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig filter --min-abundance 2 -o {output} {input}
    '''

##################################################################
## DOWNLOAD AND SKETCH PHIX
## covers:
## - phiX illumina spike in
##################################################################

rule download_phix:
    output: "outputs/sourmash_sketch/GCF_000819615.1_genomic.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz | zcat | sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund -o {output} --name GCF_000819615.1 -
    '''

##################################################################
## COMBINE TO CREATE FINAL CONTAM DATABASE
##################################################################

rule combine_and_create_contam_db:
    input: 
        db="outputs/sourmash_contam_db/gtdb-rs207.genomic-reps.dna.k{ksize}.contamsubset.zip",
        genome_sigs = expand("outputs/sourmash_sketch/{genome_acc}_cds_from_genomic.sig", genome_acc = GENOME_ACC),
        sra_sigs = expand("outputs/sourmash_sketch_filtered/{sra_acc}_sra.sig", sra_acc = SRA_ACC),
        phix_sig = "outputs/sourmash_sketch/GCF_000819615.1_genomic.sig",
    output: "outputs/sourmash_contam_db/contamdb.dna.k{ksize}.zip"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig cat -o {output} {input.db} {input.genome_sigs} {input.sra_sigs} {input.phix_sig}
    '''
