library(tidyverse)
library(uwot)

dataDir <- "../data/sdy180/resultfiles/gene_expression_result/Nanostring_norm_data_DS10_ESIDs_SDY180.587719.txt"

# ---------------------------------
# Read Data
# ---------------------------------

ns.d<-read.delim(dataDir)
head(ns.d)

d<-ns.d%>%select()
# ---------------------------------
# Convert to wide
# ---------------------------------
# Convert to wide format
d.wide <- ns.d %>% select(EXP_SAMPLE_ACC_NUM,Gene_Name,Count) %>%
  pivot_wider(names_from = Gene_Name, values_from = Count) %>%
  column_to_rownames("EXP_SAMPLE_ACC_NUM")


# ---------------------------------
# Dimensionality reduction
# ---------------------------------
# UMAP dimension reduction - repeated for each D 
	D_range <-seq(from=2,to=6)
	umap_pred<-names(d.wide)

    # Run UMAP with each D
    umap_list<-list()
    for(j in 1:length(D_range)){
      print(paste0("UMAP dimension reduction from ",paste0(length(umap_pred)),"-D to ",D_range[j],"-D"))
      set.seed(10)
      umap_out <-uwot::umap(d.wide,
                            scale=T,
                            n_threads = 10,
                            n_neighbors = 3,
                            min_dist = 0.25,
                            learning_rate = 0.5,
                            init="normlaplacian",
                            ret_model = F,
                            n_components = D_range[j],
                            verbose=T)
      
      
      umap_list[[paste0("N_dim_",D_range[j])]]<-umap_out
    }
    remove(umap_out)
    
# ---------------------------------
# Optimal BIC/max silhoutte
# ---------------------------------
# can add here to select optimal D_value for base_umap transformation
D_value<-2
suffix<-Sys.Date()
modelDir<-"/Users/croftwd/Documents/welcome_ideathon/Random_Forest_Rangers/models/"

# ---------------------------------
# Save the transform (ret_model =T)
# ---------------------------------
set.seed(10)
umap_model<-uwot::umap(d.wide,
                            scale=T,
                            n_threads = 10,
                            n_neighbors = 3,
                            min_dist = 0.25,
                            learning_rate = 0.5,
                            init="normlaplacian",
                            ret_model = T,
                            n_components = D_value,
                            verbose=T)

uwot::save_uwot(umap_model,file = paste0(modelDir,"UMAP_",D_value,"D_model_export_",suffix,".rds"))

# -------------------------------------------------------------------------------------------------------

# ---------------------------------
# Load the transform 
# ---------------------------------
 umap_model<-uwot::load_uwot(file = paste0(modelDir,"UMAP_",D_value,"D_model_export_",suffix,".rds"))

# ---------------------------------
# Clustering 
# ---------------------------------





