#########
# 01: Libraries
#########
library(tidyverse) #Dply
library(FastPG) #Clustering
library(uwot) #umap
library(FNN)  #knn
library(igraph)
library(data.table)
library(caret)
library(class)



if(!require(devtools)){
  install.packages("devtools") # If not already installed
}

devtools::install_github("JinmiaoChenLab/Rphenograph")
library(Rphenograph)
#########
# 02: Data
#########

#Training Data
sdy180 <- read.delim("./data/sdy180/resultfiles/gene_expression_result/Nanostring_norm_data_DS10_ESIDs_SDY180.587719.txt", header = TRUE, sep = "\t", dec = ".")

#Testing Data
sdy296 <- read.delim("./data/sdy296/resultfiles/gene_expression_result/Nanostring_norm_data_DS10_ESIDs_SDY296.587721.txt", header = TRUE, sep = "\t", dec = ".")

#Testing Data
sdy301 <- read.delim("./data/sdy301/resultfiles/gene_expression_result/Nanostring_norm_data_DS10_ESIDs_SDY301.587720.txt", header = TRUE, sep = "\t", dec = ".")

#########
# 03: matrix build
#########

mat_180 <- dcast(sdy180, EXP_SAMPLE_ACC_NUM ~ Gene_Name, value.var = "Count", fun.aggregate = sum)
adj_matrix <- as.matrix(mat_180[, -1])  # Remove the gene name column

mat_296 <- dcast(sdy296, EXP_SAMPLE_ACC_NUM ~ Gene_Name, value.var = "Count", fun.aggregate = sum)
mat_296 <- as.matrix(mat_296[, -1])  # Remove the gene name column

mat_301 <- dcast(sdy301, EXP_SAMPLE_ACC_NUM ~ Gene_Name, value.var = "Count", fun.aggregate = sum)
mat_301 <- as.matrix(mat_296[, -1]) 

#########
# 04: Clustering
#########

Rphenograph_out <- Rphenograph(adj_matrix, k = 5)

mat_180_umap <- umap(adj_matrix,min_dist = 0.25) 

umap_embeddings <- get_embedding(mat_180_umap)


phenograph_cluster <- factor(membership(Rphenograph_out[[2]]))
#########
# 05: UMAP (Waynes Settings)
#########

mat_180 <-  cbind(mat_180,mat_180_umap)

mat_180 <- cbind(mat_180,phenograph_cluster)

colnames(mat_180)[[62]] <- "umap_1"
colnames(mat_180)[[63]] <- "umap_2"

pdf("test_umap.pdf")
print(ggplot(mat_180,aes(x=umap_1,y=umap_2,color=factor(phenograph_cluster)))+
    geom_point())
dev.off()


#########
# 04: Training KNN model
#########

features <- mat_180[, 1:(ncol(mat_180)-3)]   # Extract the features
labels <-  mat_180[,"phenograph_cluster"]  

new_features <- mat_296

new_features <- mat_296[, intersect(colnames(mat_296), colnames(mat_180))]
features <-mat_180[, intersect(colnames(mat_180), colnames(mat_296))]


k <- 3  # Specify the number of neighbors to consider (can be adjusted as needed)

# Train the KNN model
knn_model <- knn(train = features, test =new_features, cl = labels, k = k)

#Predictions for  296 model. 
predictions <- knn(train = features, test = new_features, cl = labels, k = k)


#########
# 04: Implement 296 on umap
#########



