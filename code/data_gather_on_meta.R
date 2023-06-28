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
list.files("../data/sdy296/sdy296-dr47_tab/", full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  reduce(., union)

# find presence of column_names in each file
column_data <-
list.files(dir1, full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
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

# get intersection columns with their files
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
         size = 8,
         alpha = 0.7,
         layout.exp = 0.2)

## find path between two nodes

# link filenames with filepath
filepath_linker <-
list.files(dir1, full.names = TRUE) %>%
  data.frame(file_location = .) %>%
  mutate(filenamecsv = list.files(dir1, full.names = FALSE)) %>%
  mutate(filename = gsub("*.csv", "", filenamecsv)) 

# get edges and files
column_edges_and_files <-
combn(columns, 2, simplify = FALSE) %>%
  map(., ~get_intersect_column_plus_files(.)) %>%
  reduce(rbind) %>%
  filter(intrsct >= 1) %>%
  mutate(weight = intrsct)

# pick two random columns
nodes <- sample(columns, 2)

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
  path_files[[i]] <- read_csv(filepath_linker[filepath_linker$filename == file_graph_path[i], "file_location"], col_types = cols(.default = col_character()))
}

reduce(path_files, full_join) %>%
  print


# build it into a function
nodes <- sample(columns, 2)

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

nodes <- c("SUBJECT_ACCESSION", "EXPSAMPLE_ACCESSION")
nodes <- c("BIOSAMPLE_ACCESSION", "EXPSAMPLE_ACCESSION")
nodes <- c("ETHNICITY", "STUDY_TIME_COLLECTED")
get_data_from_columns(nodes) %>%
  print(n = 40)

full_join(path_files[[1]], path_files[[2]], na_matches = "never") %>%
  summarise_all(., ~length(levels(factor(.)))) %>%
  t()

# attempt to implement this for many column names

nodes <- sample(columns, 5)

distances(column_graph) %>%
  head