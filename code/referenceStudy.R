library(tidyverse)
library(uwot)
library(mclust)

dataFile <- "../data/sdy180/resultfiles/gene_expression_result/Nanostring_norm_data_DS10_ESIDs_SDY180.587719.txt"
expMeta<-"../data/sdy180/resultfiles/sdy180-dr47_subject_2_gene_expression_result.txt"

# ---------------------------------
# Read Data
# ---------------------------------

ns.d<-read.delim(dataFile)
head(ns.d)

ns.meta<-read.delim(expMeta)
# ---------------------------------
# Convert to wide
# ---------------------------------
# Convert to wide format
d.wide <- ns.d %>% select(EXP_SAMPLE_ACC_NUM,Gene_Name,Count) %>%
  pivot_wider(names_from = Gene_Name, values_from = Count) %>%
  column_to_rownames("EXP_SAMPLE_ACC_NUM")%>%select(grep("^NEG",value=TRUE,invert=TRUE,names(.)))


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
D_value<-2 # We are setting manually to 2
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
saveRDS(GMM_model,paste0(modelDir,"GMM_k_",K_value,"_D_",D_value,"_model_",suffix,".rds"))
    

# ---------------------------------
# Basic UMAP plots
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

embed<-embed%>%rownames_to_column("Expsample.Accession")%>%left_join(.,ns.meta,by="Expsample.Accession")

ggArm<-ggplot(embed,aes_string(x="X1",y="X2", colour="ARM.Accession"))+
        geom_point(size=3,alpha=0.7)+
        theme_bw() +
        labs(color="Study Arm")+
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              legend.position="top")+
        xlab("UMAP1") + ylab("UMAP2")

ggGender<-ggplot(embed,aes_string(x="X1",y="X2", colour="Gender"))+
        geom_point(size=3,alpha=0.7)+
        theme_bw() +
        labs(color="Gender")+
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              legend.position="top")+
        xlab("UMAP1") + ylab("UMAP2")

ggAge<-ggplot(embed,aes_string(x="X1",y="X2", colour="Subject.Age"))+
        geom_point(size=3,alpha=0.7)+
        theme_bw() +
        labs(color="Age")+
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              legend.position="top")+
        xlab("UMAP1") + ylab("UMAP2")

ggTime<-ggplot(embed,aes_string(x="X1",y="X2", colour="Planned.Visit.Name"))+
        geom_point(size=3,alpha=0.7)+
        theme_bw() +
        labs(color="Timepoint")+
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              legend.position="top")+
        xlab("UMAP1") + ylab("UMAP2")


ggGender + ggAge + ggArm + ggTime
ggsave("../results/SDY180_ns_umap.png",width=8,height=2.2)






