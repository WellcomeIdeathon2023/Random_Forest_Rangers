library(tidyverse)
library(igraph)
library(GGally)
library(RColorBrewer)
library(knitr)
library(magrittr)

dir1 <- "../data/sdy296/sdy296-dr47_tab/"

# get file names in directory
file_names <-
list.files(dir1) %>%
  gsub("*.csv", "", .) %>%
  as.list()

# get column names (includes ...1)
column_names <-
list.files(dir1, full.names = TRUE) %>%
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

# construct a graph to analyse file relationships

# remove rownumber column
link <-
column_data %>%
  filter(column != "...1") %>%
  filter(present) 

# find number of connecting files
get_intersect_file <- function(x) {
  a1 <- pull(filter(link, file == x[[1]]), column)
  b1 <- pull(filter(link, file == x[[2]]), column)

  data.frame(a = x[[1]],
             b = x[[2]],
             intrsct = length(intersect(a1, b1))
             ) %>% return
}

# find number of connecting files and return the files
get_intersect_column_plus_files <- function(x) {
  a1 <- pull(filter(link, column == x[[1]]), file)
  b1 <- pull(filter(link, column == x[[2]]), file)

  tibble(a = x[[1]],
             b = x[[2]],
             intrsct = length(intersect(a1, b1)),
             files = list(intersect(a1, b1))
             ) %>% return
}

# get file names
files <-
link %>%
  pull(file) %>%
  unique

# compare each file to see intersect of presence in ACCESSION columns
file_edges <-
combn(files, 2, simplify = FALSE) %>%
  map(., ~get_intersect_file(.)) %>%
  reduce(rbind) %>%
  filter(intrsct >= 1) %>%
  mutate(weight = intrsct)

# create graph
file_graph <-
file_edges %>%
  graph_from_data_frame(directed = FALSE)


# cluster the graph
file_clusters <-
cluster_louvain(file_graph, resolution = 1.2)

# plot

generate_hex_colours <- function(x){
  brewer.pal(length(unique(x)), "Dark2")[as.numeric(factor(x))]
}

file_graph %>%
  ggnet2(label = TRUE,
         label.size = 3,
         label.alpha = 0.6,
         color = file_clusters$membership,
         size = 20,
         alpha = 0.3,
         layout.exp = 0.2)


# link filenames with filepath
filepath_linker <-
list.files(dir1, full.names = TRUE) %>%
  data.frame(file_location = .) %>%
  mutate(filenamecsv = list.files(dir1, full.names = FALSE)) %>%
  mutate(filename = gsub("*.csv", "", filenamecsv)) 

## find path between two files and return dataframes
file_nodes <- sample(files, 2)

get_data_from_files <- function(nodes){
  path1 <-
  shortest_paths(file_graph, from = file_nodes[[1]],
                 to = file_nodes[[2]],
                 output = "vpath") %>%
    extract2("vpath") %>% extract2(1) %>% names(.)

  path_files <- list()
  for(i in 1:length(path1)){
    path_files[[i]] <- read_csv(filepath_linker[filepath_linker$filename == path1[i], "file_location"], col_types = cols(.default = col_character())) %>%
      select(., -...1)
  }

  reduce(path_files, full_join) %>%
    return
}

# find which keys to join on (using ACCESSION)
joining_keys <-
intersect(colnames(path_files[[2]]), colnames(path_files[[3]])) %>%
  .[grepl("ACCESSION", .)]

# other columns which have the same name need to be removed or else clashes occur
non_joining_keys <-
intersect(colnames(path_files[[2]]), colnames(path_files[[3]])) %>%
  .[!grepl("ACCESSION", .)]

full_join(path_files[[2]], path_files[[3]])

a <- path_files[[2]][,!colnames(path_files[[2]]) %in% non_joining_keys]

b <- path_files[[3]][,!colnames(path_files[[3]]) %in% non_joining_keys]

full_join(a, b,
          na_matches = "never",
          by = setNames(joining_keys, joining_keys)
          )

path_files[[2]]$NAME
path_files[[3]]$NAME

# this function will take a vector of two column names only
get_data_from_columns <- function(nodes){
  # find shortest path between columns
  path1 <-
  shortest_paths(column_graph, from = nodes[[1]],
                 to = nodes[[2]],
                 output = "vpath") %>%
    extract2("vpath") %>% extract2(1) %>% names(.)

  # loop over path to find files in order
  file_graph_path <- vector()
  for(i in 1:(length(path1)-1)){
    # filter needs to work in both directions as entries not duplicated in a/b
    file_graph_path <-
    filter(column_edges_and_files,
           a %in% path1[c(i, i+1)],
           b %in% path1[c(i, i+1)]
           ) %>%
    pull(files) %>% extract2(1) %>%
    append(file_graph_path, .)
  }

  path_files <- list()
  for(i in 1:length(file_graph_path)){
    path_files[[i]] <- read_csv(filepath_linker[filepath_linker$filename == file_graph_path[i], "file_location"], col_types = cols(.default = col_character())) %>%
      select(., -...1)
  }

  reduce(path_files, full_join) %>%
    return
}

# pick two random columns
nodes <- sample(columns, 2)

get_data_from_columns(nodes) %>%
  select(nodes, everything()) %>%
  print
  write.csv("../results/merge_test.csv")

# so do all files have an _ACCESSION

# get column names (includes ...1)
column_names_all <-
list.files(dir1, full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  reduce(., union) %>%
  .[!grepl("...1", .)]

# find presence of column_names in each file
list.files(dir1, full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  map(., ~grepl("ACCESSION", .)) %>%
  map(., ~sum(.)) %>%
  reduce(., c)
# all files contain an _ACCESSION column

# get all columns presence in files
column_data_all <-
list.files(dir1, full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  map2(., file_names, ~data.frame(column = column_names, file = .y, present = column_names %in% .x)) %>%
  reduce(rbind) %>%
  filter(column != "...1") 

nodes <- sample(column_names_all, 2)

# find which files have these data present

files_of_interest <-
column_data_all %>%
  filter(column %in% nodes) %>%
  filter(present) %>%
  pull(file)




# attempt to implement this for many column names

nodes <- sample(columns, 5)

distances(column_graph) %>%
  head
