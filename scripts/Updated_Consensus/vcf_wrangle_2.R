# Wrangle variant calling data - version 2 - using merged VCF file

# Load Packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(tibble)

VCF <- read_tsv("/home/jwestaway/pk_pipeline/Initial_Pk_subset/outputs/variant_calling/consensus/updated/merged_variants_only.tsv", col_names = TRUE) %>% 
  rename("CHROM" = "#CHROM") %>% 
  mutate_at(1:5, as.factor) %>% 
  unite(variant_name, CHROM, POS, ID, REF, ALT, sep = "_") %>% 
  select(-c(2:5)) %>% 
  mutate_at(c(2:ncol(.)), ~str_remove(., ":.*")) %>% 
  mutate_at(c(2:ncol(.)), ~ifelse(. %like% "1/1", "1", .)) %>% # homozygous
  mutate_at(c(2:ncol(.)), ~ifelse(. %like% "1/0" | . %like% "0/1" , "0.5", .)) %>%  # heterozygous
  mutate_at(c(2:ncol(.)), ~ifelse(. != "1" & . != "0.5" , "0", .)) %>% 
  mutate(variant_name = str_remove(variant_name, "_\\.")) 

ordered_VCF <- VCF %>% 
  filter(grepl("ordered", variant_name)) %>% 
  select(1, starts_with("2:")) %>% 
  pivot_longer(!variant_name, names_to = "sample", values_to = "variant_1") %>% 
  mutate(sample = str_remove(sample, "2:")) %>% 
  left_join(
    VCF %>% 
      filter(grepl("ordered", variant_name)) %>% 
      select(1, !starts_with("2:")) %>% 
      pivot_longer(!variant_name, names_to = "sample", values_to = "variant_2")
  ) %>% 
  mutate_at(c(3:4), as.numeric) %>% 
  mutate(consensus_variant = ifelse((variant_1 == variant_2) & (variant_1 + variant_2 >= 1), 1, 0)) %>%
  filter(consensus_variant == 1) %>% 
  select(1,2,5) %>% 
  pivot_wider(names_from = sample, values_from = consensus_variant) %>%  
  select(1) %>%  
  separate(variant_name, into = c("X1", "X2", "X3", "X4", "POS", "REF", "ALT"), sep = "_") %>%  
  unite(CHROM, starts_with("X")) %>% 
  add_column(ID = ".") %>% 
  relocate(CHROM, POS, ID)  

MIT_VCF <- VCF %>% 
  filter(grepl("MIT", variant_name)) %>% 
  select(1, starts_with("2:")) %>% 
  pivot_longer(!variant_name, names_to = "sample", values_to = "variant_1") %>% 
  mutate(sample = str_remove(sample, "2:")) %>% 
  left_join(
    VCF %>% 
      filter(grepl("MIT", variant_name)) %>% 
      select(1, !starts_with("2:")) %>% 
      pivot_longer(!variant_name, names_to = "sample", values_to = "variant_2")
  ) %>% 
  mutate_at(c(3:4), as.numeric) %>% 
  mutate(consensus_variant = ifelse((variant_1 == variant_2) & (variant_1 + variant_2 >= 1), 1, 0)) %>%
  filter(consensus_variant == 1) %>% 
  select(1,2,5) %>% 
  pivot_wider(names_from = sample, values_from = consensus_variant) %>%  
  select(1) %>%  
  separate(variant_name, into = c("X1", "X2", "X3", "POS", "REF", "ALT"), sep = "_") %>%  
  unite(CHROM, starts_with("X")) %>% 
  add_column(ID = ".") %>% 
  relocate(CHROM, POS, ID)  

API_VCF <- VCF %>% 
  filter(grepl("API", variant_name)) %>% 
  select(1, starts_with("2:")) %>% 
  pivot_longer(!variant_name, names_to = "sample", values_to = "variant_1") %>% 
  mutate(sample = str_remove(sample, "2:")) %>% 
  left_join(
    VCF %>% 
      filter(grepl("API", variant_name)) %>% 
      select(1, !starts_with("2:")) %>% 
      pivot_longer(!variant_name, names_to = "sample", values_to = "variant_2")
  ) %>% 
  mutate_at(c(3:4), as.numeric) %>% 
  mutate(consensus_variant = ifelse((variant_1 == variant_2) & (variant_1 + variant_2 >= 1), 1, 0)) %>%
  filter(consensus_variant == 1) %>% 
  select(1,2,5) %>% 
  pivot_wider(names_from = sample, values_from = consensus_variant) %>%  
  select(1) %>%  
  separate(variant_name, into = c("X1", "X2", "X3", "X4", "X5", "POS", "REF", "ALT"), sep = "_") %>%  
  unite(CHROM, starts_with("X")) %>% 
  add_column(ID = ".") %>% 
  relocate(CHROM, POS, ID)  

ordered_VCF %>% 
  rbind(MIT_VCF, API_VCF) %>% 
  write_tsv("/home/jwestaway/pk_pipeline/Initial_Pk_subset/outputs/variant_calling/consensus/updated/vcf_variant_names.tsv")

