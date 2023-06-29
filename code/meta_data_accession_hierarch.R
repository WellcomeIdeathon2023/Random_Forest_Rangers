library(tidyverse)
library(igraph)
library(GGally)
library(RColorBrewer)
library(knitr)

dir1 <- "../data/sdy296/sdy296-dr47_tab/"

# get file names in directory
file_names <-
list.files(dir1) %>%
  gsub("*.csv", "", .) %>%
  as.list()

# link filenames with filepath
filepath_linker <-
list.files(dir1, full.names = TRUE) %>%
  data.frame(file_location = .) %>%
  mutate(filenamecsv = list.files(dir1, full.names = FALSE)) %>%
  mutate(filename = gsub("*.csv", "", filenamecsv)) 

# get column names (includes ...1)
column_names <-
list.files("../data/sdy296/sdy296-dr47_tab/", full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  reduce(., union) %>%
  .[grepl("ACCESSION", .)]

# find presence of column_names in each file
column_data <-
list.files(dir1, full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  map(., ~(.[grepl("ACCESSION", .)])) %>%
  map2(., file_names, ~data.frame(column = column_names, file = .y, present = column_names %in% .x)) %>%
  reduce(rbind) 

# construct a graph to analyse column relationships

# remove rownumber column
link <-
column_data %>%
  filter(column != "...1") %>%
  filter(present) 

# find number of intersecting files per column
get_intersect_column <- function(x) {
  a1 <- pull(filter(link, column == x[[1]]), file)
  b1 <- pull(filter(link, column == x[[2]]), file)

  data.frame(a = x[[1]],
             b = x[[2]],
             intrsct = length(intersect(a1, b1))
             ) %>% return
}

get_intersect_column_plus_files <- function(x) {
  a1 <- pull(filter(link, column == x[[1]]), file)
  b1 <- pull(filter(link, column == x[[2]]), file)

  tibble(a = x[[1]],
             b = x[[2]],
             intrsct = length(intersect(a1, b1)),
             files = list(intersect(a1, b1))
             ) %>% return
}

# get column names without ...1
columns <-
link %>%
  pull(column) %>%
  unique

# compare each column to see intersect of presence in files
column_edges <-
combn(columns, 2, simplify = FALSE) %>%
  map(., ~get_intersect_column(.)) %>%
  reduce(rbind) %>%
  filter(intrsct >= 1) %>%
  mutate(weight = intrsct)

# analyse each edge to determine direction

column_edges_and_files <-
combn(columns, 2, simplify = FALSE) %>%
  map(., ~get_intersect_column_plus_files(.)) %>%
  reduce(rbind) %>%
  filter(intrsct >= 1) %>%
  mutate(weight = intrsct)

link1 <- column_edges_and_files[1,]

b <- link1$files[[1]][[1]]

column_mapping <-
read_csv(filepath_linker[filepath_linker$filename == b, "file_location"]) %>%
  select(link1$a, link1$b) %>%
  distinct %>%
  map(., ~dup_both(.)) %>%
  map(., ~sum(.)) %>%
  map(., ~ifelse(. > 0, "MANY", "ONE")) 


dup_both <- function(x) duplicated(x) | duplicated(x, fromLast = TRUE)

dup_both(c(1, 2, 2, 1, 2, 2, 3, 1, 3, 4))


# create graph
column_graph <-
column_edges %>%
  graph_from_data_frame(directed = FALSE)


# cluster the graph
column_clusters <-
cluster_louvain(column_graph, resolution = 1.2)

# plot
generate_hex_colours <- function(x){
  brewer.pal(length(unique(x)), "Dark2")[as.numeric(factor(x))]
}

column_graph %>%
  ggnet2(label = TRUE,
         label.size = 3,
         label.alpha = 0.6,
         color = column_clusters$membership,
         size = 20,
         alpha = 0.3,
         layout.exp = 0.2)

ggsave("../results/sdy296_meta_accession_weighted_clustered_graph.png", height = 8, width = 8)
