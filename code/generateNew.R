library(tidyverse)

modelDir<-"../models/"
study<-c("SDY296","SDY301")
dataFiles <- c("../data/sdy296/resultfiles/gene_expression_result/Nanostring_norm_data_DS10_ESIDs_SDY296.587721.txt",
  "../data/sdy301/resultfiles/gene_expression_result/Nanostring_norm_data_DS10_ESIDs_SDY301.587720.txt")

metaFiles<-c("../data/sdy180/resultfiles/sdy180-dr47_subject_2_gene_expression_result.txt",
  "../data/sdy296/resultfiles/sdy296-dr47_subject_2_gene_expression_result.txt",
  "../data/sdy301/resultfiles/sdy301-dr47_subject_2_gene_expression_result.txt")

cluster_mod<-list.files(modelDir,pattern="GMM*")
umap_mod<-list.files(modelDir,pattern="UMAP*")

# ---------------------------------
# Read Data
# ---------------------------------
ns.d<-lapply(dataFiles,read.delim)
names(ns.d)<-study

meta.list<-lapply(metaFiles,read.delim)
meta.d<-bind_rows(meta.list[[1]],meta.list[[2]],meta.list[[3]])
meta.d<-meta.d%>%select(.,c("Subject.Accession","Gender","Subject.Age","ARM.Accession","Expsample.Accession"))

# ---------------------------------
# Convert to wide
# ---------------------------------
include<-c("Gender","Subject.Age") # Additional data to include 

d.wide<-list()
for(i in study){
  d.wide[[i]] <- ns.d[[i]] %>% select(EXP_SAMPLE_ACC_NUM,Gene_Name,Count) %>%
  pivot_wider(names_from = Gene_Name, values_from = Count) %>%
  select(grep("^NEG",value=TRUE,invert=TRUE,names(.)))%>%
  select(grep("^POS",value=TRUE,invert=TRUE,names(.)))%>%
  rename(Expsample.Accession=EXP_SAMPLE_ACC_NUM)%>%
  left_join(.,meta.d[,c("Expsample.Accession",include)],by="Expsample.Accession")%>%
  column_to_rownames("Expsample.Accession")
}

# ---------------------------------
# New patient example randomly generated from ranges in SDY296
# For demonstrating utility in predicting new data 
# ---------------------------------
# Example data frame
df <- d.wide[[1]][,1:46]
df$Subject.Age<-d.wide[[1]]$Subject.Age
# Function to generate random numbers within the range of a column
generate_random_row <- function(col) {
  min_val <- min(col)
  max_val <- max(col)
  random_val <- runif(1, min_val, max_val)
  return(random_val)
}

# Add a new row with random numbers within the range of each column
new_sample <- apply(df, 2, generate_random_row)
data.frame(new_sample)
new_sample$Gender=sample(c('Male', 'Female'),1)
new_sample$Subject.Age<-round(new_sample$Subject.Age)
new_sample<-data.frame(new_sample[names(d.wide[[1]])])
rownames(new_sample)<-"new_sample"

# Write to file
write.csv(new_sample,"../data/randomNew.csv")




