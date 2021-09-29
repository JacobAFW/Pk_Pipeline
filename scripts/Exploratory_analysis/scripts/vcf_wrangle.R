# Wrangle variant calling data

# Load Packages
library(tidyverse)

# Read in and wrangle data to get full summary data

## Need to create names for the columns as the first row only has 7 columns, and thus R assumes all rows only have 7 columns and we end up losing a significant amount of data
Names <- c("X1") # X1 is the first columns
for (i in 1:100){
  Names <- append(Names, print(paste0("Sample_", i)))
}

Summary_Data <- read_tsv("C:/Users/Jacob/Desktop/Menzies/Knowlesi/Pipeline/Pk_Pipeline/data/vcf/merged.tsv", col_names =  Names) %>% 
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
  mutate_at(c("DP", "GQ", "MQ"), as.numeric)

  write_tsv(Summary_Data, "vcf_full_summary.csv")

# Read in and wrangle data to summaries of each - the n() will actually give the number of samples with variants
Summary_Data %>% 
  unite(Variant, Contig, Base, ID, Ref, Alt, sep = "-") %>% # create a unique variant code (Variant) based on these variables
  group_by(V_Call_Tool) %>% 
  dplyr::summarise(Samples_with_variants = n(), 
                   GQ = mean(GQ),
                   DP = mean(DP),
                   MQ = mean(MQ)) %>% 
  write_tsv("vcf_summary.csv")
  

# Get counts of Variants 
Summary_Data %>% 
  unite(Variant, Contig, Base, ID, Ref, Alt, sep = "-") %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(V_Call_Tool) %>% 
  dplyr::summarise(Variants = length(unique(Variant))) %>% 
  write_tsv("vcf_variant_counts_in_tool.csv")


# Get variants that are found within both tools
Summary_Data %>% 
  unite(Variant, Contig, Base, ID, Ref, Alt, sep = "-") %>%  # create a unique variant code (Variant) based on these variables
  select(Variant, V_Call_Tool) %>% 
  filter(V_Call_Tool == "BCFTOOLS") %>% 
  mutate_if(is.character, as.factor) %>% 
  unique() %>% 
  inner_join(
    (Summary_Data %>% 
      unite(Variant, Contig, Base, ID, Ref, Alt, sep = "-") %>%  
      select(Variant, V_Call_Tool) %>% 
      filter(V_Call_Tool == "GATK") %>% 
      mutate_if(is.character, as.factor) %>% 
      unique()), by = "Variant") %>% 
  select("Variant") %>% 
  separate(Variant, sep = "-", c("CHROM", "POS", "ID", "REF", "ALT")) %>% 
  write_tsv("vcf_variant_names.csv")


# Get the FMT scores for the vairants found in both tools
Summary_Data %>% 
  unite(Variant, Contig, Base, ID, Ref, Alt, sep = "-") %>% # create a unique variant code (Variant) based on these variables
  filter(V_Call_Tool == "BCFTOOLS") %>% 
  select(!c(Sample, PL, V_Call_Tool)) %>% 
  group_by(Variant) %>% 
  summarise_all(mean) %>%  # mean value for each varint
  inner_join(
    (Summary_Data %>% 
       unite(Variant, Contig, Base, ID, Ref, Alt, sep = "-") %>% 
       filter(V_Call_Tool == "GATK") %>% 
       select(!c(Sample, PL, V_Call_Tool)) %>% 
       group_by(Variant) %>% 
       summarise_all(mean)), 
    by = "Variant") %>% # inner join by variants to only select those the are found by both tools
  mutate(DP = (DP.x + DP.y)/2, # calculating the means for the two tools
         GQ = (GQ.x + GQ.y)/2,
         MQ = (MQ.x + MQ.y)/2) %>% 
  select(1, 8:10) %>% 
  write_tsv("vcf_variants_and_metrics.csv")








