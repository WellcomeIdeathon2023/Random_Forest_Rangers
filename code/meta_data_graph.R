library(tidyverse)
library(igraph)
library(GGally)

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

# remove rownumber column
link <-
column_data %>%
  filter(column != "...1") %>%
  filter(present) 

# find number of intersecting column names
get_intersect_file <- function(x) {
  a1 <- pull(filter(link, file == x[[1]]), column)
  b1 <- pull(filter(link, file == x[[2]]), column)

  data.frame(a = x[[1]],
             b = x[[2]],
             intrsct = length(intersect(a1, b1))
             ) %>% return
}

# get files
files <-
column_data %>%
  pull(file) %>%
  unique

# compare each file with each other and find intersection of column names
file_edges <-
combn(files, 2, simplify = FALSE) %>%
  map(., ~get_intersect_file(.)) %>%
  reduce(rbind) %>%
  filter(intrsct >= 1)

# create non-directional unweighted graph
file_graph <-
file_edges %>%
  graph_from_data_frame(directed = FALSE)

# plot graph using fructermanreingold
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
  filter(intrsct >= 1)

# create graph
column_graph <-
column_edges %>%
  graph_from_data_frame(directed = FALSE)

# plot
column_graph %>%
  ggnet2(label = TRUE,
         label.size = 3,
         label.alpha = 0.6,
         size = 8,
         color = "#c99cc0",
         alpha = 0.7,
         layout.exp = 0.2)

ggsave("../results/sdy296_variable_relationship_graph.png", height = 8, width = 8)
