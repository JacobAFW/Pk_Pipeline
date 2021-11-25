# Wrangle variant calling data

# Load Packages
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(forcats)

Names <- c("X1") # X1 is the first columns
for (i in 1:100){
  Names <- append(Names, print(paste0("Sample_", i)))
}

read_tsv("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/GATK_query.tsv", col_names =  Names) %>% 
  select(X1) %>% 
  inner_join(read_tsv("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/samtools_query.tsv", col_names =  Names) %>% 
               select(X1)) %>% 
  separate(X1, sep =" ", c("Contig", "Base", "ID", "Ref", "Alt")) %>% 
  write_tsv("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/vcf_variant_names.tsv")