library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(forcats)
library(ggplot2)

read_tsv("runs_of_homozygosity_head.tsv", col_names = F) %>% 
  ggplot(mapping = aes(x=X4, y=X6)) + 
  geom_point()
  
ggsave("runs_of_homo.png")