# This script filters signatures in a sourmash database to those that have been detected as common kit contaminants in doi 10.1016/j.tim.2018.11.003 Table1.
# Species Klebsiella aerogenes is also added as this is an organisms being used by Arcadians.
# To accomplish this subsetting, the genuses of the contaminants are read in and the GTDB metadata is filtered to species in the data base that have those genuses.
# Only representatives for that species are retained to keep the database small (e.g. we only need one e. coli genome, not 5000).

library(readr)
library(dplyr)
library(tidyr)

bac120 <- read_tsv(snakemake@input[['gtdb_metadata']]) %>%
  filter(gtdb_representative == TRUE) %>%
  separate(gtdb_taxonomy, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep  = ";")

contam <- read_csv(snakemake@input[['contams']]) %>%
  mutate(genus = paste0("g__", genus)) %>%
  mutate(genus = gsub("g__Enhydrobacter", "g__Moraxella_A", genus))  %>% # g__Enhydrobacter is g__Moraxella_A in GTDB
  mutate(genus = gsub("g__kingella", "g__Kingella", genus)) # fix capitalization so it matches

# check transfer criteria
# contam$genus[!contam$genus %in% bac120$genus]
# TM7 isn't really a taxonomic lineage, so left out

# subset the gtdb metadata to potential contaminant organisms
bac_contam1 <- bac120 %>%
  filter(genus %in% contam$genus)

bac_contam2 <- bac120 %>%
  filter(species %in% c("s__Klebsiella aerogenes"))

bac_contam <- bind_rows(bac_contam1, bac_contam2)

# remove GB and RF prefixes from accessions 
bac_contam <- bac_contam %>%
  mutate(accession = gsub("^GB_", "", accession),
         accession = gsub("^RS_", "", accession))

# read in sig describe and subset to genomes in bac_contam data frame
sig_describe <- read_csv(snakemake@input[['sig_describe']]) %>%
  mutate(genome_acc = gsub(" .*", "", name)) %>%
  filter(genome_acc %in% bac_contam$accession) %>%
  select(-genome_acc)

# write out the outputs
write_csv(sig_describe, snakemake@output[["picklist"]])
