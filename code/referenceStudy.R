library(tidyverse)
library(uwot)
library(mclust)

dataFile <- "../data/sdy180/resultfiles/gene_expression_result/Nanostring_norm_data_DS10_ESIDs_SDY180.587719.txt"

# ---------------------------------
# Read Data
# ---------------------------------

ns.d<-read.delim(dataFile)
head(ns.d)

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

# ---------------------------------
# GMM model-based clustering on reference dataset UMAP transform
# ---------------------------------
# GMM model & save
suffix <- Sys.Date()
# Choose K, the number of components in the model
K_value = 3 # could iterate and optimise best k
# Choose constraint for GMM covariance (note, unconstrained covariance matrices will dramatically increase run time)
mod_type = "VVV" 
    
# Base umap
dat_umap <- umap_model$embedding
    
# GMM mclust
#.Random.seed <- seed_save
set.seed(10) 
GMM_model <- Mclust(umap_model$embedding,
                        G = K_value,
                        modelNames = mod_type,
                        initialization = list("hcpairs"))

# Clusters
GMM_model$classification

# Saving GMM model
saveRDS(GMM_model,paste0(modelDir,"GMM_k_",K_value,"D_",D_value,"_model_",suffix,".rds"))
    

# ---------------------------------
# Basic UMAP plot
# ---------------------------------
embed<-data.frame(umap_model$embedding,stringsAsFactors = F)

gg<-ggplot(embed,aes_string(x="X1",y="X2"))+
      geom_point(size=1,alpha=0.7)+
      labs(x="UMAP1", y="UMAP2") +
      theme_bw() +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())#+
      #scale_color_manual(values=c("#00BFC4","#F8766D"))
    gg






