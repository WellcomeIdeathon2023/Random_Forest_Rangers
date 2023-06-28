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

columns <-
list.files("../data/sdy296/sdy296-dr47_tab/", full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  reduce(., union)


file_names <-
list.files("../data/sdy296/sdy296-dr47_tab/") %>%
  gsub("*.csv", "", .) %>%
  as.list()

column_data <-
list.files("../data/sdy296/sdy296-dr47_tab/", full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  map2(., file_names, ~data.frame(column = columns, file = .y, present = columns %in% .x)) %>%
  reduce(rbind) 


column_order <-
column_data %>%
  group_by(column) %>%
  summarise(n = sum(present)) %>%
  arrange(desc(n)) %>%
  pull(column)

name_order <-
column_data %>%
  group_by(file) %>%
  summarise(n = sum(present)) %>%
  arrange(desc(n)) %>%
  pull(file)

column_data %>%
  group_by(file) %>%
  mutate(column_n = sum(present)) %>%
  ungroup %>%
  group_by(column) %>%
  mutate(file_n = sum(present)) %>%
  ungroup %>%
  arrange(desc(column_n), desc(file_n)) %>%
  select(-column_n, -file_n) %>%
  pivot_wider(names_from = c("file"), values_from = "present", names_sort = FALSE) %>%
  write.csv("../results/sdy296_meta_data.csv")

column_data %>%
  group_by(file) %>%
  mutate(column_n = sum(present)) %>%
  ungroup %>%
  group_by(column) %>%
  mutate(file_n = sum(present)) %>%
  ungroup %>%
  arrange(desc(column_n), desc(file_n)) %>%
  select(-column_n, -file_n) %>%
  pivot_wider(names_from = c("file"), values_from = "present", names_sort = FALSE) %>%
  print(n = 50)

# attempt to join all colums
# list.files("../data/sdy296/sdy296-dr47_tab/", full.names = TRUE) %>%
#   map(., ~read_csv(., col_types = cols(.default = col_character()))) %>%
#   map(., ~select(., -...1)) %>%
#   reduce(full_join)
# this fails because it is not joining in order

# construct a graph to analyse file relationships

library(igraph)
library(GGally)

link <-
column_data %>%
  filter(column != "...1") %>%
  filter(present) 

get_intersect_file <- function(x) {
  a1 <- pull(filter(link, file == x[[1]]), column)
  b1 <- pull(filter(link, file == x[[2]]), column)

  data.frame(a = x[[1]],
             b = x[[2]],
             intrsct = length(intersect(a1, b1))
             ) %>% return
}

files <-
column_data %>%
  pull(file) %>%
  unique

file_edges <-
combn(files, 2, simplify = FALSE) %>%
  map(., ~get_intersect_file(.)) %>%
  reduce(rbind) %>%
  filter(intrsct >= 1)

file_graph <-
file_edges %>%
  graph_from_data_frame(directed = FALSE)

file_graph %>%
  ggnet2(label = TRUE,
         label.size = 3,
         label.alpha = 0.6,
         size = 25,
         color = "#9cc0c9",
         alpha = 0.7,
         layout.exp = 0.2)

ggsave("../results/sdy296_file_relationship_graph.png", width = 6, height = 6)


# construct a graph to analyse column relationships

link <-
column_data %>%
  filter(column != "...1") %>%
  filter(present) 

get_intersect_column <- function(x) {
  a1 <- pull(filter(link, column == x[[1]]), file)
  b1 <- pull(filter(link, column == x[[2]]), file)

  data.frame(a = x[[1]],
             b = x[[2]],
             intrsct = length(intersect(a1, b1))
             ) %>% return
}

columns <-
link %>%
  pull(column) %>%
  unique

column_edges <-
combn(columns, 2, simplify = FALSE) %>%
  map(., ~get_intersect_column(.)) %>%
  reduce(rbind) %>%
  filter(intrsct >= 1)

column_graph <-
column_edges %>%
  graph_from_data_frame(directed = FALSE)

column_graph %>%
  ggnet2(label = TRUE,
         label.size = 3,
         label.alpha = 0.6,
         size = 8,
         color = "#c99cc0",
         alpha = 0.7,
         layout.exp = 0.2)

ggsave("../results/sdy296_variable_relationship_graph.png", height = 8, width = 8)
