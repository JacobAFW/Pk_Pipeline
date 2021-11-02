# Load Packages
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(forcats)

#########################################################################################################################################################

# Define function for variant totals

variant_count_total <- function(INDEL_TSV, SNP_TSV, FILTER){

Names <- c("X1") # X1 is the first columns
for (i in 1:100){
  Names <- append(Names, print(paste0("Sample_", i)))
}
# Need to create names for the columns as the first row only has 7 columns, and thus R assumes all rows only have 7 columns and we end up losing a significant amount of data
Variants <- read_tsv(INDEL_TSV, col_names =  Names) %>% 
  separate(X1, sep =" ", c("Contig", "Base", "ID", "Ref", "Alt")) %>%  
  pivot_longer(cols = !c(Contig, Base, ID, Ref, Alt)) %>%  
  select(-name) %>% 
  na.omit() %>% 
  mutate(V_Call_Tool = ifelse(grepl("2:", value), "BCFTOOLS", "GATK")) %>% 
  separate(value, sep = " ", c("Sample", "DP", "GQ", "MQ", "PL")) %>% 
  mutate(DP = str_remove(DP, "DP=")) %>% 
  mutate(GQ = str_remove(GQ, "GQ=")) %>% 
  mutate(MQ = str_remove(MQ, "MQ=")) %>% 
  mutate(PL = str_remove(PL, "PL=")) %>% 
  mutate_at(c("DP", "GQ", "MQ"), as.numeric) %>% 
  add_column(Variant = "Indel") %>% 
  rbind(
   read_tsv(SNP_TSV, col_names =  Names) %>% 
     separate(X1, sep =" ", c("Contig", "Base", "ID", "Ref", "Alt")) %>% 
     pivot_longer(cols = !c(Contig, Base, ID, Ref, Alt)) %>%  
     select(-name) %>% 
     na.omit() %>% 
     mutate(V_Call_Tool = ifelse(grepl("2:", value), "BCFTOOLS", "GATK")) %>% 
     separate(value, sep = " ", c("Sample", "DP", "GQ", "MQ", "PL")) %>% 
     mutate(DP = str_remove(DP, "DP=")) %>% 
     mutate(GQ = str_remove(GQ, "GQ=")) %>% 
     mutate(MQ = str_remove(MQ, "MQ=")) %>% 
     mutate(PL = str_remove(PL, "PL=")) %>% 
     mutate_at(c("DP", "GQ", "MQ"), as.numeric) %>% 
     add_column(Variant = "SNP")  
   ) %>% 
  separate(Alt, c("ALT1", "ALT2", "ALT3", "ALT4", "ALT5", "ALT6"), sep = ",") %>% 
  pivot_longer(5:10, names_to = "ALT_N", values_to = "Alt") %>% 
  select(-ALT_N) %>% 
  na.omit(Alt)

Variants_2 <- Variants %>% # up to here produces a df that lists every variant-sample combination & below summarises this to give us counts
  unite(Variant_ID, Contig, Base, ID, Ref, Alt, sep = "-") %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(Variant) %>% 
  dplyr::summarise(Variant_Count = length(unique(Variant_ID)), 
                   GQ = mean(GQ),
                   DP = mean(DP),
                   MQ = mean(MQ)) 

Variants_2 %>% 
  add_column(Filter = FILTER) %>% 
  select(-c(GQ, DP, MQ)) %>% 
  pivot_wider(c(Filter, Variant_Count), names_from = Variant, values_from = Variant_Count) %>% 
  left_join(
    Variants_2 %>% 
  add_column(Filter = FILTER) %>% 
    select(Filter, GQ, DP, MQ) %>% 
    group_by(Filter) %>% 
    summarise_all(mean)) %>% 
  mutate(Total = Indel + SNP) %>% 
  relocate(Filter, Total, SNP, Indel)
}

#########################################################################################################################################################

# Define function for tool-specific counts

variant_count_tool_spec <- function(INDEL_TSV, SNP_TSV, FILTER, TOOL){

Names <- c("X1") # X1 is the first columns
for (i in 1:100){
  Names <- append(Names, print(paste0("Sample_", i)))
}
# Need to create names for the columns as the first row only has 7 columns, and thus R assumes all rows only have 7 columns and we end up losing a significant amount of data
Variants <- read_tsv(INDEL_TSV, col_names =  Names) %>% 
  separate(X1, sep =" ", c("Contig", "Base", "ID", "Ref", "Alt")) %>% 
  pivot_longer(cols = !c(Contig, Base, ID, Ref, Alt)) %>%  
  select(-name) %>% 
  na.omit() %>% 
  separate(value, sep = " ", c("Sample", "DP", "GQ", "MQ", "PL")) %>% 
  mutate(DP = str_remove(DP, "DP=")) %>% 
  mutate(GQ = str_remove(GQ, "GQ=")) %>% 
  mutate(MQ = str_remove(MQ, "MQ=")) %>% 
  mutate(PL = str_remove(PL, "PL=")) %>% 
  mutate_at(c("DP", "GQ", "MQ"), as.numeric) %>% 
  add_column(Variant = "Indel") %>% 
  rbind(
   read_tsv(SNP_TSV, col_names =  Names) %>% 
     separate(X1, sep =" ", c("Contig", "Base", "ID", "Ref", "Alt"))  %>% 
     pivot_longer(cols = !c(Contig, Base, ID, Ref, Alt)) %>%  
     select(-name) %>% 
     na.omit() %>% 
     separate(value, sep = " ", c("Sample", "DP", "GQ", "MQ", "PL")) %>% 
     mutate(DP = str_remove(DP, "DP=")) %>% 
     mutate(GQ = str_remove(GQ, "GQ=")) %>% 
     mutate(MQ = str_remove(MQ, "MQ=")) %>% 
     mutate(PL = str_remove(PL, "PL=")) %>% 
     mutate_at(c("DP", "GQ", "MQ"), as.numeric) %>% 
     add_column(Variant = "SNP")  
   ) %>% 
  separate(Alt, c("ALT1", "ALT2", "ALT3", "ALT4", "ALT5", "ALT6"), sep = ",") %>% 
  pivot_longer(5:10, names_to = "ALT_N", values_to = "Alt") %>% 
  select(-ALT_N) %>% 
  na.omit(Alt)

Variants_2 <- Variants %>% # up to here produces a df that lists every variant-sample combination & below summarises this to give us counts
  unite(Variant_ID, Contig, Base, ID, Ref, Alt, sep = "-") %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(Variant) %>% 
  dplyr::summarise(Variant_Count = length(unique(Variant_ID)), 
                   GQ = mean(GQ),
                   DP = mean(DP),
                   MQ = mean(MQ)) 

Variants_2 %>% 
  add_column(Filter = FILTER) %>% 
  select(-c(GQ, DP, MQ)) %>% 
  pivot_wider(c(Filter, Variant_Count), names_from = Variant, values_from = Variant_Count) %>% 
  left_join(
    Variants_2 %>% 
  add_column(Filter = FILTER) %>% 
    select(Filter, GQ, DP, MQ) %>% 
    group_by(Filter) %>% 
    summarise_all(mean)) %>% 
  mutate(Total = Indel + SNP) %>% 
  relocate(Filter, Total, SNP, Indel)  %>% 
  rename_at(vars(SNP), funs(paste0(TOOL, "_SNP"))) %>% 
  rename_at(vars(Indel), funs(paste0(TOOL, "_Indel"))) %>% 
  select(1,3:4)

}

#########################################################################################################################################################

# Use functions with left_join to combine 'orgiginal data'

# ORIGINAL
ORIGINAL <- variant_count_total("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/original/indels.tsv", 
                    "/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/original/SNPs.tsv", 
                    "ORIGINAL") %>% 
  left_join(variant_count_tool_spec("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/original/indels_GATK.tsv", 
                        "/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/original/SNPs_GATK.tsv", 
                        "ORIGINAL", "GATK")) %>% 
  left_join(variant_count_tool_spec("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/original/indels_bcftools.tsv", 
                        "/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/original/SNPs_bcftools.tsv", 
                        "ORIGINAL", "bcftools")) 


#########################################################################################################################################################

# Alter Functions for filters

#########################################################################################################################################################

# Define function for variant totals

variant_count_total <- function(INDEL_TSV, SNP_TSV, FILTER){

Names <- c("X1") # X1 is the first columns
for (i in 1:100){
  Names <- append(Names, print(paste0("Sample_", i)))
}
# Need to create names for the columns as the first row only has 7 columns, and thus R assumes all rows only have 7 columns and we end up losing a significant amount of data
Variants <- read_tsv(INDEL_TSV, col_names =  Names) %>% 
  separate(X1, sep =" ", c("Contig", "Base", "ID", "Ref", "Alt", "FILTER")) %>% 
  select(-FILTER) %>% 
  pivot_longer(cols = !c(Contig, Base, ID, Ref, Alt)) %>%  
  select(-name) %>% 
  na.omit() %>% 
  mutate(V_Call_Tool = ifelse(grepl("2:", value), "BCFTOOLS", "GATK")) %>% 
  separate(value, sep = " ", c("Sample", "DP", "GQ", "MQ", "PL")) %>% 
  mutate(DP = str_remove(DP, "DP=")) %>% 
  mutate(GQ = str_remove(GQ, "GQ=")) %>% 
  mutate(MQ = str_remove(MQ, "MQ=")) %>% 
  mutate(PL = str_remove(PL, "PL=")) %>% 
  mutate_at(c("DP", "GQ", "MQ"), as.numeric) %>% 
  add_column(Variant = "Indel") %>% 
  rbind(
   read_tsv(SNP_TSV, col_names =  Names) %>% 
     separate(X1, sep =" ", c("Contig", "Base", "ID", "Ref", "Alt", "FILTER")) %>% 
     select(-FILTER) %>% 
     pivot_longer(cols = !c(Contig, Base, ID, Ref, Alt)) %>%  
     select(-name) %>% 
     na.omit() %>% 
     mutate(V_Call_Tool = ifelse(grepl("2:", value), "BCFTOOLS", "GATK")) %>% 
     separate(value, sep = " ", c("Sample", "DP", "GQ", "MQ", "PL")) %>% 
     mutate(DP = str_remove(DP, "DP=")) %>% 
     mutate(GQ = str_remove(GQ, "GQ=")) %>% 
     mutate(MQ = str_remove(MQ, "MQ=")) %>% 
     mutate(PL = str_remove(PL, "PL=")) %>% 
     mutate_at(c("DP", "GQ", "MQ"), as.numeric) %>% 
     add_column(Variant = "SNP")  
   ) %>% 
  separate(Alt, c("ALT1", "ALT2", "ALT3", "ALT4", "ALT5", "ALT6"), sep = ",") %>% 
  pivot_longer(5:10, names_to = "ALT_N", values_to = "Alt") %>% 
  select(-ALT_N) %>% 
  na.omit(Alt)

Variants_2 <- Variants %>% # up to here produces a df that lists every variant-sample combination & below summarises this to give us counts
  unite(Variant_ID, Contig, Base, ID, Ref, Alt, sep = "-") %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(Variant) %>% 
  dplyr::summarise(Variant_Count = length(unique(Variant_ID)), 
                   GQ = mean(GQ),
                   DP = mean(DP),
                   MQ = mean(MQ)) 

Variants_2 %>% 
  add_column(Filter = FILTER) %>% 
  select(-c(GQ, DP, MQ)) %>% 
  pivot_wider(c(Filter, Variant_Count), names_from = Variant, values_from = Variant_Count) %>% 
  left_join(
    Variants_2 %>% 
  add_column(Filter = FILTER) %>% 
    select(Filter, GQ, DP, MQ) %>% 
    group_by(Filter) %>% 
    summarise_all(mean)) %>% 
  mutate(Total = Indel + SNP) %>% 
  relocate(Filter, Total, SNP, Indel)
}

#########################################################################################################################################################

# Define function for tool-specific counts

variant_count_tool_spec <- function(INDEL_TSV, SNP_TSV, FILTER, TOOL){

Names <- c("X1") # X1 is the first columns
for (i in 1:100){
  Names <- append(Names, print(paste0("Sample_", i)))
}
# Need to create names for the columns as the first row only has 7 columns, and thus R assumes all rows only have 7 columns and we end up losing a significant amount of data
Variants <- read_tsv(INDEL_TSV, col_names =  Names) %>% 
  separate(X1, sep =" ", c("Contig", "Base", "ID", "Ref", "Alt", "FILTER")) %>% 
  select(-FILTER) %>% 
  pivot_longer(cols = !c(Contig, Base, ID, Ref, Alt)) %>%  
  select(-name) %>% 
  na.omit() %>% 
  separate(value, sep = " ", c("Sample", "DP", "GQ", "MQ", "PL")) %>% 
  mutate(DP = str_remove(DP, "DP=")) %>% 
  mutate(GQ = str_remove(GQ, "GQ=")) %>% 
  mutate(MQ = str_remove(MQ, "MQ=")) %>% 
  mutate(PL = str_remove(PL, "PL=")) %>% 
  mutate_at(c("DP", "GQ", "MQ"), as.numeric) %>% 
  add_column(Variant = "Indel") %>% 
  rbind(
   read_tsv(SNP_TSV, col_names =  Names) %>% 
     separate(X1, sep =" ", c("Contig", "Base", "ID", "Ref", "Alt", "FILTER")) %>% 
     select(-FILTER) %>% 
     pivot_longer(cols = !c(Contig, Base, ID, Ref, Alt)) %>%  
     select(-name) %>% 
     na.omit() %>% 
     separate(value, sep = " ", c("Sample", "DP", "GQ", "MQ", "PL")) %>% 
     mutate(DP = str_remove(DP, "DP=")) %>% 
     mutate(GQ = str_remove(GQ, "GQ=")) %>% 
     mutate(MQ = str_remove(MQ, "MQ=")) %>% 
     mutate(PL = str_remove(PL, "PL=")) %>% 
     mutate_at(c("DP", "GQ", "MQ"), as.numeric) %>% 
     add_column(Variant = "SNP")  
   ) %>% 
  separate(Alt, c("ALT1", "ALT2", "ALT3", "ALT4", "ALT5", "ALT6"), sep = ",") %>% 
  pivot_longer(5:10, names_to = "ALT_N", values_to = "Alt") %>% 
  select(-ALT_N) %>% 
  na.omit(Alt)

Variants_2 <- Variants %>% # up to here produces a df that lists every variant-sample combination & below summarises this to give us counts
  unite(Variant_ID, Contig, Base, ID, Ref, Alt, sep = "-") %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(Variant) %>% 
  dplyr::summarise(Variant_Count = length(unique(Variant_ID)), 
                   GQ = mean(GQ),
                   DP = mean(DP),
                   MQ = mean(MQ)) 

Variants_2 %>% 
  add_column(Filter = FILTER) %>% 
  select(-c(GQ, DP, MQ)) %>% 
  pivot_wider(c(Filter, Variant_Count), names_from = Variant, values_from = Variant_Count) %>% 
  left_join(
    Variants_2 %>% 
  add_column(Filter = FILTER) %>% 
    select(Filter, GQ, DP, MQ) %>% 
    group_by(Filter) %>% 
    summarise_all(mean)) %>% 
  mutate(Total = Indel + SNP) %>% 
  relocate(Filter, Total, SNP, Indel)  %>% 
  rename_at(vars(SNP), funs(paste0(TOOL, "_SNP"))) %>% 
  rename_at(vars(Indel), funs(paste0(TOOL, "_Indel"))) %>% 
  select(1,3:4)

}

#########################################################################################################################################################

# Use functions with left_join to combine filters

# GATK_BP
GATK_BP <- variant_count_total("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/GATK_BP/indels.tsv", 
                    "/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/GATK_BP/SNPs.tsv", 
                    "GATK_BP") %>% 
  left_join(variant_count_tool_spec("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/GATK_BP/indels_GATK.tsv", 
                        "/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/GATK_BP/SNPs_GATK.tsv", 
                        "GATK_BP", "GATK")) %>% 
  left_join(variant_count_tool_spec("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/GATK_BP/indels_bcftools.tsv", 
                        "/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/GATK_BP/SNPs_bcftools.tsv", 
                        "GATK_BP", "bcftools")) 

# VIVAX
VIVAX <- variant_count_total("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/Vivax/indels.tsv", 
                    "/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/Vivax/SNPs.tsv", 
                    "VIVAX") %>% 
  left_join(variant_count_tool_spec("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/Vivax/indels_GATK.tsv", 
                        "/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/Vivax/SNPs_GATK.tsv", 
                        "VIVAX", "GATK")) %>% 
  left_join(variant_count_tool_spec("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/Vivax/indels_bcftools.tsv", 
                        "/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/Vivax/SNPs_bcftools.tsv", 
                        "VIVAX", "bcftools")) 

#########################################################################################################################################################

# Use rbind to create a summary files

ORIGINAL %>%
    rbind(GATK_BP) %>%
    rbind(VIVAX) %>%
    write_tsv("/home/jwestaway/pk_pipeline/ZB_100/outputs/vcf_consensus/filtering_tests/Filtering_Summary.tsv")
