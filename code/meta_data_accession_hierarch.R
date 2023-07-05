library(tidyverse)
library(magrittr)
library(igraph)
library(GGally)
library(RColorBrewer)
library(knitr)
library(cowplot)

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

# find presence of column_names all in each file
column_names_all <-
list.files("../data/sdy296/sdy296-dr47_tab/", full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  reduce(., union) 

# find presence of column_names in each file
column_data <-
list.files(dir1, full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  map(., ~(.[grepl("ACCESSION", .)])) %>%
  map2(., file_names, ~data.frame(column = column_names, file = .y, present = column_names %in% .x)) %>%
  reduce(rbind) 

# presence of all column_names in each file
column_data_all <-
list.files(dir1, full.names = TRUE) %>%
  map(., ~read_csv(.)) %>%
  map(., ~colnames(.)) %>%
  map2(., file_names, ~data.frame(column = column_names_all, file = .y, present = column_names_all %in% .x)) %>%
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

dup_both <- function(x) duplicated(x) | duplicated(x, fromLast = TRUE)

column_edges_and_files <-
combn(columns, 2, simplify = FALSE) %>%
  map(., ~get_intersect_column_plus_files(.)) %>%
  reduce(rbind) %>%
  filter(intrsct >= 1) %>%
  mutate(weight = intrsct)

column_edges_and_files <- mutate(column_edges_and_files, a_mapping = NA, b_mapping = NA)

for (j in 1:nrow(column_edges_and_files)){
  link1 <- column_edges_and_files[j,]

  column_mappings <- list()
  for(i in 1:length(link1$files[[1]])){

    b <- link1$files[[1]][[i]]

    column_mapping <-
    read_csv(filepath_linker[filepath_linker$filename == b, "file_location"]) %>%
      select(a_mapping = link1$a, b_mapping = link1$b) %>%
      distinct() %>%
      map(., ~dup_both(.)) %>%
      map(., ~sum(.)) %>%
      map(., ~ifelse(. > 0, "DUP", "UNIQUE")) 

    column_mappings[[i]] <- data.frame(column_mapping)
  }

  #  amalgamate analysis into one line
  column_mapping_amal <-
  reduce(column_mappings, rbind) %>%
    distinct 

  # check that we have one mapping for all files
  if(nrow(column_mapping_amal) == 1){
    column_edges_and_files[j, "a_mapping"] <- column_mapping_amal$a_mapping
    column_edges_and_files[j, "b_mapping"] <- column_mapping_amal$b_mapping
  } else {
  # if files are discrepant, then "MANY" takes priority
  # i.e. one-to-many should overrule one-to-one.
    column_mapping_amal <- 
    column_mapping_amal %>%
      map(., ~ifelse(sum(. == "DUP") > 0, "DUP", "UNIQUE")) %>%
      data.frame(.)
    message("mappings not unique")
    column_edges_and_files[j, "a_mapping"] <- column_mapping_amal$a_mapping
    column_edges_and_files[j, "b_mapping"] <- column_mapping_amal$b_mapping
  }
}

column_edge_map_type <-
column_edges_and_files %>%
  mutate(map_type = ifelse(a_mapping != b_mapping, "one-to-many",
                           ifelse(a_mapping == "UNIQUE", "one-to-one", "many-to-many"))) %>%
  select(-files) 

column_edge_map_type %>%
  ggplot(aes(x = map_type))+
  geom_bar()

# swap around values in vertices to provide direction for igraph
column_edge_map_type_directed <- column_edge_map_type
for (i in 1:nrow(column_edge_map_type)){
  if(column_edge_map_type[i,"map_type"] == "one-to-many" & column_edge_map_type[i, "a_mapping"] == "DUP"){
    b <- column_edge_map_type[i, "a"]
    a <- column_edge_map_type[i, "b"]
    b_mapping <- column_edge_map_type[i, "a_mapping"]
    a_mapping <- column_edge_map_type[i, "b_mapping"]

    column_edge_map_type_directed[i, "a"] <- a
    column_edge_map_type_directed[i, "b"] <- b
    column_edge_map_type_directed[i, "a_mapping"] <- a_mapping
    column_edge_map_type_directed[i, "b_mapping"] <- b_mapping
  }
  if(column_edge_map_type[i, "map_type"] %in% c("one-to-one", "many-to-many")){
    b <- column_edge_map_type[i, "a"]
    a <- column_edge_map_type[i, "b"]
    b_mapping <- column_edge_map_type[i, "a_mapping"]
    a_mapping <- column_edge_map_type[i, "b_mapping"]

    column_edge_map_type_directed[i, "a"] <- a
    column_edge_map_type_directed[i, "b"] <- b
    column_edge_map_type_directed[i, "a_mapping"] <- a_mapping
    column_edge_map_type_directed[i, "b_mapping"] <- b_mapping

    # duplicate
    column_edge_map_type_directed <- rbind(column_edge_map_type_directed, column_edge_map_type[i,])
  }
}

# plot
generate_hex_colours <- function(x){
  brewer.pal(length(unique(x)), "Dark2")[as.numeric(factor(x))]
}

# create graph
column_graph <-
column_edge_map_type_directed %>%
  mutate(color = generate_hex_colours(map_type)) %>%
  graph_from_data_frame(directed = TRUE)

column_graph %>%
  ggnet2(label = TRUE,
         mode = "fruchtermanreingold",
         label.size = 3,
         label.alpha = 0.6,
         edge.color = "color",
         size = 20,
         alpha = 0.3,
         arrow.size = 5,
         arrow.gap = 0.02,
         layout.exp = 0.2)

ggsave("../results/sdy296_meta_data_directed_graph.png", height = 8, width = 8)

# algorithm to score the nodes
edge_path_find <- column_edge_map_type_directed
for (i in 1:nrow(column_edge_map_type_directed)){
  if(column_edge_map_type_directed[i,"map_type"] == "one-to-many" & column_edge_map_type_directed[i, "a_mapping"] == "UNIQUE"){
    b <- column_edge_map_type_directed[i, "a"]
    a <- column_edge_map_type_directed[i, "b"]
    b_mapping <- column_edge_map_type_directed[i, "a_mapping"]
    a_mapping <- column_edge_map_type_directed[i, "b_mapping"]

    edge_path_find[i, "a"] <- a
    edge_path_find[i, "b"] <- b
    edge_path_find[i, "a_mapping"] <- a_mapping
    edge_path_find[i, "b_mapping"] <- b_mapping
  }
}

# need direction to be ONE -> MANY
# also don't want ONE - ONE, or MANY - MANY

edge_path_directed <-
edge_path_find %>%
  filter(map_type == "one-to-many")

get_next_node <- function(vert1){
edge_path_directed %>%
  filter(a == vert1) %>%
  sample_n(., 1) %>%
  .$b %>%
  return
}

# run random walks
iterations <- 300

{
  vertices_score <- data.frame(vert = unique(c(edge_path_directed$a, edge_path_directed$b)))

  vert1 <- sample(vertices_score$vert, 1)

  # start at 5
  vertices_score[vertices_score$vert == vert1, "score"] <- 5
  adjusted <- 1
  run <- TRUE

  while (run) {
    iterations = iterations - 1
    if (iterations < 1) {run = FALSE}

    # check if vert1 is na, if so, assign it a number
    vert1_score <- vertices_score[vertices_score$vert == vert1, "score"]
    if (is.na(vert1_score)) {
      vert1_score <- 5}
    vertices_score[vertices_score$vert == vert1, "score"] <- vert1_score 

    # check if there is a vert2 that it points to
    if (nrow(filter(edge_path_directed, a == vert1)) > 0){

      # get the vert2
      vert2 <- get_next_node(vert1)

      # get vert2 score
      vert2_score <- vertices_score[vertices_score$vert == vert2, "score"]

      # compare score with vert1 and adjust if needed
      if (is.na(vert2_score)) {
        vert2_score <- vert1_score + 1 
        adjusted <- adjusted + 1
      } else if (vert2_score <= vert1_score) {
        vert2_score <- vert1_score + 1
        adjusted <- adjusted + 1
      } else if (vert2_score >  vert1_score) {
      }

      # write out new vert2 score
      vertices_score[vertices_score$vert == vert2, "score"] <- vert2_score 

      # get next vert2
      vert1 <- vert2

    } else {


      # if there is no vert2 that it points to, generate a new vert1 for next cycle
      if (sum(is.na(vertices_score$score)) > 0) {
        vert1 <-
          vertices_score %>%
          filter(is.na(score)) %>%
          sample_n(., 1) %>%
          .$vert
      } else {
        vert1 <- vertices_score %>% sample_n(., 1) %>% .$vert
      }
    }

    print(adjusted)
  }
}

vert_scor_joined <- left_join(data.frame(vert = names(V(column_graph))), vertices_score)

#visual check on graph scores
column_graph %>%
  ggnet2(label = vert_scor_joined$score,
         mode = "fruchtermanreingold",
         label.size = 3,
         label.alpha = 0.6,
         edge.color = "color",
         size = 20,
         alpha = 0.3,
         arrow.size = 5,
         arrow.gap = 0.02,
         layout.exp = 0.2)

ggsave("../results/sdy_296_meta_istrahler_graph.png")

vertices_score %>%
  arrange(desc(score))

# get scores
scores <-
vertices_score %>%
  pull(score) %>%
  unique

file_graphs <- list()
file_edges_plus_cols <- list()
for (i in 1:length(scores)) {

  # get column names
  columns_names_subset <-
  vertices_score %>%
    filter(score >= scores[i]) %>%
    pull(vert)

  # find presence of column_names in each file
  column_data_subset <-
  list.files(dir1, full.names = TRUE) %>%
    map(., ~read_csv(.)) %>%
    map(., ~colnames(.)) %>%
    map(., ~.[. %in% columns_names_subset]) %>%
    map2(., file_names, ~data.frame(column = column_names, file = .y, present = column_names %in% .x)) %>%
    reduce(rbind) 

  # construct a graph to analyse column relationships

  # remove rownumber column
  link <-
  column_data_subset %>%
    filter(column != "...1") %>%
    filter(present) 

  # find number of intersecting columns per file
  get_intersect_file <- function(x) {
    a1 <- pull(filter(link, file == x[[1]]), column)
    b1 <- pull(filter(link, file == x[[2]]), column)

    data.frame(a = x[[1]],
               b = x[[2]],
               intrsct = length(intersect(a1, b1))
               ) %>% return
  }

  get_intersect_file_plus_column <- function(x) {
    a1 <- pull(filter(link, file == x[[1]]), column)
    b1 <- pull(filter(link, file == x[[2]]), column)

    tibble(a = x[[1]],
               b = x[[2]],
               intrsct = length(intersect(a1, b1)),
               columns = list(intersect(a1, b1))
               ) %>% return
  }

  # get file names
  files <-
  link %>%
    pull(file) %>%
    unique

  # compare each column to see intersect of presence in files
  file_edges <-
  combn(files, 2, simplify = FALSE) %>%
    map(., ~get_intersect_file(.)) %>%
    reduce(rbind) %>%
    filter(intrsct >= 1) %>%
    mutate(weight = intrsct)

  file_edges_plus_cols[[i]] <-
  combn(files, 2, simplify = FALSE) %>%
    map(., ~get_intersect_file_plus_column(.)) %>%
    reduce(rbind) %>%
    filter(intrsct >= 1) %>%
    mutate(weight = intrsct)

  # create graph
  file_graphs[[i]] <-
  file_edges %>%
    graph_from_data_frame(directed = FALSE)
}

plots <-
file_graphs %>%
  map(., ~ggnet2(., label = TRUE,
         mode = "fruchtermanreingold",
         label.size = 3,
         label.alpha = 0.6,
         size = 20,
         alpha = 0.3,
         arrow.size = 5,
         arrow.gap = 0.02,
         layout.exp = 0.2)
  )

file_edges_plus_cols[[2]]

cowplot::plot_grid(plotlist = plots)

# choose data needed
column_names_all
columns_needed <- c("EXPERIMENT_ACCESSION", "ETHNICITY")

# find the files they are present in 
files_needed <-
column_data_all %>%
  filter(column %in% columns_needed) %>%
  filter(present) %>%
  pull(file)

# join files by highest scoring graph
files_to_join <- files_needed[1:2]

# find both files in graph

for (i in length(file_edges_plus_cols):1){

  graph_i <- file_edges_plus_cols[[i]]

  present_in_graph <-
  graph_i %>%
    filter(a %in% files_to_join & b %in% files_to_join) %>%
    nrow() %>% is_greater_than(0)

  if(present_in_graph){
    break
  }
}

if(present_in_graph){
  print("present in graph")
} else {
  print("not found in graph")
}

# join graphs

if(present_in_graph){
  join_by_columns <-
  graph_i %>%
    filter(a %in% files_to_join & b %in% files_to_join) %>%
    pull(columns) %>%
    extract2(1)

    extract2(columns)

    print


# create data creating algorithm
# steps:
# 1. choose data you want
# 2. locate files that contain data
# 3. create multiple graph of files based on vert score
# 4. start with highest scoring graph
# 5. can you connect files using this graph?
# 6. if not proceed down graph
# 7. when you can, then you have wondering tradesman problem
# 8. join by best route
