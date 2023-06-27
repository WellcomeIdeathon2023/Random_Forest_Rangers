library(tidyverse)

# SDY296

# number of genes in RNA seq
read_tsv("../data/sdy296/resultfiles/rna_sequencing_result/SDY296_EXP13760_RNA_seq.703270.tsv") %>%
  pivot_longer(-c("ENSEMBL", "SYMBOL", "TYPE")) %>%
  pull(ENSEMBL) %>% unique %>% length

individuals <- list()

# number of individuals in nanostring
individuals[["ge"]] <-
read_tsv("../data/sdy296/resultfiles/sdy296-dr47_subject_2_gene_expression_result.txt") %>%
  pull(`Subject Accession`) %>%
  unique

# number of individuals in rna seq
individuals[["rna"]] <-
read_tsv("../data/sdy296/resultfiles/sdy296-dr47_subject_2_rna_sequencing_result.txt") %>%
  pull(`Subject Accession`) %>%
  unique

# number of individuals in rna seq
individuals[["hai"]] <-
read_csv("../data/sdy296/resultfiles/hai_result.csv") %>%
  pull(SUBJECT_ACCESSION) %>%
  unique

# number of individuals in neut_ab
individuals[["neut_ab"]] <-
read_csv("../data/sdy296/resultfiles/neut_ab_titer_result.csv") %>%
  pull(SUBJECT_ACCESSION) %>%
  unique

reduce(individuals, union) # 37 total patients
reduce(individuals, intersect) # 13 patients with data in all 4

# number of HAI serology
read_csv("../data/sdy296/resultfiles/hai_result.csv") %>%
  pull(SUBJECT_ACCESSION) %>%
  unique %>% length # n = 37

read_csv("../data/sdy296/resultfiles/neut_ab_titer_result.csv") %>%
  pull(SUBJECT_ACCESSION) %>%
  unique %>% length # n = 37

