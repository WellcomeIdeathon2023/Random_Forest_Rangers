# Uniform Manifold Approximation and Projection -assisted clustering workflow:
# --------------------------------------------------------------------------------------------------------------------
# referenceStudy.R
# 	- reads in Nanostring expression data
# 	- increased dimensionality is challenging for clustering algorithms so compression and noise reduction helps
# 	- capture nonlinearities with UMAP dimensionality reduction to 2D (can use mixed numerical and categorical data so include age/demographics etc)
# 	- identify subgroups that predict protection
# 	- returns the model allowing new data to be added to an existing embedding
# 	- can optimise D selection with scree plot and BIC metrics

library(tidyverse)
library(uwot)
library(mclust)
library(ComplexHeatmap)
library(patchwork)
library(RColorBrewer)
library(circlize)
library(ggpubr)

dataFile <- "../data/sdy180/resultfiles/gene_expression_result/Nanostring_norm_data_DS10_ESIDs_SDY180.587719.txt"
expMeta<-"../data/sdy180/resultfiles/sdy180-dr47_subject_2_gene_expression_result.txt"

# ---------------------------------
# Read Data
# ---------------------------------

ns.d<-read.delim(dataFile)
head(ns.d)

ns.meta<-read.delim(expMeta)
# ---------------------------------
# Convert ns data to wide
# ---------------------------------
# Convert to wide format
d.wide <- ns.d %>% select(EXP_SAMPLE_ACC_NUM,Gene_Name,Count) %>%
  pivot_wider(names_from = Gene_Name, values_from = Count) %>%
  column_to_rownames("EXP_SAMPLE_ACC_NUM")%>%select(grep("^NEG",value=TRUE,invert=TRUE,names(.)))

# ---------------------------------
# Add additional variables to include
# ---------------------------------
include<-c("Gender","Subject.Age") # Additional data to include 

d.wide <- ns.d %>% select(EXP_SAMPLE_ACC_NUM,Gene_Name,Count) %>%
  pivot_wider(names_from = Gene_Name, values_from = Count) %>%
  select(grep("^NEG",value=TRUE,invert=TRUE,names(.)))%>% # remove negative control probes
  select(grep("^POS",value=TRUE,invert=TRUE,names(.)))%>% # remove positive control probes
  rename(Expsample.Accession=EXP_SAMPLE_ACC_NUM)%>%
  left_join(.,ns.meta[,c("Expsample.Accession",include)],by="Expsample.Accession")%>%
  column_to_rownames("Expsample.Accession")


# ---------------------------------
# Dimensionality reduction
# ---------------------------------
# UMAP dimension reduction - repeated for each D 
# Currently we manually fix D but could iterate to find optimum
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

modelDir<-"../models/"

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

if(!file.exists(paste0(modelDir,"UMAP_",D_value,"D_model_export.rds"))){
	uwot::save_uwot(umap_model,file = paste0(modelDir,"UMAP_",D_value,"D_model_export.rds"))
}

# ---------------------------------
# GMM model-based clustering on reference dataset UMAP transform
# ---------------------------------
# GMM model & save
# Choose K, the number of components in the model
K_value = 3 # could iterate and optimise best k
# Choose constraint for GMM covariance (note, unconstrained covariance matrices will dramatically increase run time)
mod_type = "VVV" # multivariate mixture ellipsoidal, varying vol shape and orientation
      
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
saveRDS(GMM_model,paste0(modelDir,"GMM_k_",K_value,"_D_",D_value,"_model_.rds"))
    

# ---------------------------------
# Basic UMAP plots
# ---------------------------------


# Colour Palette ---------------------------------------------------
pal<-list()
pal$Group<-c("#32CD32","#F37FB8","#409388","#CE8BAE","#B23648",
                  "#ADD8E6","#D46E7E","#7E486B","#79AA7A","#FFEC2B",
                  "#8D5B96","#E41A1C","#00B4F0","#3A85A8","#488846",
                  "#BD6253","#46A169","#EB7AA9","#C4625D","#D8B62E",
                  "#d6c624","#77777C","#4F6FA1","#E1712E","#A65628",
                  "#B392A3","#E984B9","#F2E631","#999999")
pal$Study<-c(SDY180="brown",SDY296="purple", SDY301="grey")
pal$Arm<-c(ARM773="#E984B9",ARM776="#F2E631",ARM779="#999999")
pal$Time<-c("Day_-7"="#C4625D","Day_0"="#D8B62E")
pal$Gender<-c(Female="red",Male="blue")
names(pal$Group) <- seq(from=1,to=length(pal$Group),by=1)
# ------------------------------------------------------------------

embed<-data.frame(umap_model$embedding,stringsAsFactors = F)

gg<-ggplot(embed,aes_string(x="X1",y="X2"))+
      geom_point(size=5,alpha=0.7)+
      labs(x="UMAP1", y="UMAP2") +
      theme_bw() +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank())
      #scale_color_manual(values=c("#00BFC4","#F8766D"))
    gg

    ggsave("../results/umap_dimension_reduction.png")

embed<-embed%>%rownames_to_column("Expsample.Accession")%>%left_join(.,ns.meta,by="Expsample.Accession")%>%
  mutate(Group=factor(GMM_model$classification))%>%mutate(Study.Time.Collected=paste0("Day_",Study.Time.Collected))


ggArm<-ggplot(embed,aes_string(x="X1",y="X2", colour="ARM.Accession"))+
        geom_point(size=3,alpha=0.7)+
        theme_bw() +
        labs(color="Study Arm")+
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              legend.position="top")+
        xlab("UMAP1") + ylab("UMAP2")+
        scale_color_manual(values=pal$Arm)

ggGender<-ggplot(embed,aes_string(x="X1",y="X2", colour="Gender"))+
        geom_point(size=3,alpha=0.7)+
        theme_bw() +
        labs(color="Gender")+
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              legend.position="top")+
        xlab("UMAP1") + ylab("UMAP2")+
        scale_color_manual(values=pal$Gender)

ggAge<-ggplot(embed,aes_string(x="X1",y="X2", colour="Subject.Age"))+
        geom_point(size=3,alpha=0.7)+
        theme_bw() +
        labs(color="Age")+
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              legend.position="top")+
        xlab("UMAP1") + ylab("UMAP2")

ggTime<-ggplot(embed,aes_string(x="X1",y="X2", colour="Study.Time.Collected"))+
        geom_point(size=3,alpha=0.7)+
        theme_bw() +
        labs(color="Timepoint")+
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              legend.position="top")+
        xlab("UMAP1") + ylab("UMAP2")+
        scale_color_manual(values=pal$Time)


ggGender + ggAge + ggArm + ggTime
ggsave("../results/SDY180_ns_umap.png",width=8,height=7)

# ---------------------------------
# Cluster characterisations
# ---------------------------------

# need to set annotation colour palette
annotation_colors <- list(#Subject.Accession=palette.top$highLevelTypePal,
                          ARM.Accession=pal$Arm,
                          Study.Accession=pal$Study,
                          Study.Time.Collected=pal$Time,
                          Group=pal$Group)

d<-d.wide%>%select(-include)
row_ha = rowAnnotation(count = anno_boxplot(d,which="row",size = unit(0.5, "mm"), width = unit(1, "cm"),box_width = 0.3))

ha = HeatmapAnnotation(df = embed[,c("Subject.Accession","ARM.Accession","Study.Accession","Study.Time.Collected","Group")],
                       show_annotation_name = TRUE,
                       col = annotation_colors,
                       simple_anno_size = unit(0.3, "cm"),
                       show_legend = c(FALSE,TRUE,TRUE,TRUE,TRUE))

# Expression data
my_data <- t(as.matrix(d))
my_data<-t(scale(t(my_data)))
# Heatmap
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
png(file = "../results/SDY180_ns_heat.png",
    width = 500,
    height = 500,)
Heatmap(
  my_data,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_order = NULL,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  bottom_annotation = NULL,
  row_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(title = ""),
  top_annotation = ha,
  left_annotation = row_ha
)
dev.off()

# --------------------------------------------------------------------------------------------------------------------

# umap_GMM_predict.R
# 	- Unsupervised clustering uses probabilistic model that assumes data is generated from a mixture of Gaussian distributions
# 	- can optimise K (number of mixture components) by silhouette width

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
# Load the transform and clustering model 
# ---------------------------------
umap_model<-uwot::load_uwot(file = paste0(modelDir,umap_mod))
gmm_model<-readRDS(paste0(modelDir,cluster_mod))

# ---------------------------------
# Check conistent variables are present in new data
# ---------------------------------

# ---------------------------------
# Transform and cluster new observations
# ---------------------------------
UMAP_t<-list()
pred<-list()
set.seed(10)
for(i in study){
  UMAP_t[[i]] <- uwot::umap_transform(d.wide[[i]],umap_model)
  # Make cluster p/Users/croftwd/Documents/welcome_ideathon/Random_Forest_Rangers/models/redictions using GMM model built on reference data
  pred[[i]] <- predict.Mclust(newdata = UMAP_t[[i]], # This transformation 
                               object = gmm_model) # Fit model 
names(pred[[i]]$classification)<-rownames(UMAP_t[[i]])
}

# ---------------------------------
# Plot reference and new data clusters
# ---------------------------------

embed<-data.frame(umap_model$embedding,stringsAsFactors = F,
  "Expsample.Accession"=rownames(umap_model$embedding),
  "Group"=factor(gmm_model$classification))
embed$SDY="SDY180"

embed<-inner_join(embed,meta.list[[1]][,c("Expsample.Accession","Subject.Accession","Gender","Subject.Age")],by="Expsample.Accession")%>%column_to_rownames("Expsample.Accession")

colour_vec<-c("#32CD32","#F37FB8","#409388","#CE8BAE","#B23648",
                  "#ADD8E6","#D46E7E","#7E486B","#79AA7A","#FFEC2B",
                  "#8D5B96","#E41A1C","#00B4F0","#3A85A8","#488846",
                  "#BD6253","#46A169","#EB7AA9","#C4625D","#D8B62E",
                  "#d6c624","#77777C","#4F6FA1","#E1712E","#A65628",
                  "#B392A3","#E984B9","#F2E631","#999999")
names(colour_vec) <- seq(from=1,to=length(colour_vec),by=1)

ggRef<-ggplot(embed,aes_string(x="X1",y="X2", colour="Group"))+
      geom_point(size=3,alpha=0.7)+
      theme_bw() +
      labs(color="GMM cluster")+
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank())+
      xlab("UMAP1") + ylab("UMAP2")+
      scale_color_manual(values=colour_vec) + ggtitle("SDY180")
    ggRef

# New data UMAPS with cluster overlay
ggNew<-list()
embedNew<-list()
newStudyMeta<-meta.list[2:3]
names(newStudyMeta)<-study
for(i in study){
    e<-as.data.frame(UMAP_t[[i]])%>%mutate(Expsample.Accession=rownames(.))%>%
    mutate(Group=factor(pred[[i]]$classification))%>%rename(,X1=V1)%>%rename(,X2=V2)%>%mutate(SDY=i)
  
   e<-inner_join(e,newStudyMeta[[i]][,c("Expsample.Accession","Subject.Accession","Gender","Subject.Age")],by="Expsample.Accession")%>%column_to_rownames("Expsample.Accession")

   embedNew[[i]]<-e

   ggNew[[i]]<-ggplot(embedNew[[i]],aes_string(x="X1",y="X2", colour="Group"))+
        geom_point(size=3,alpha=0.7)+
        theme_bw() +
        labs(color="GMM cluster")+
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              legend.position="none")+
        xlab("UMAP1") + ylab("UMAP2")+
        scale_color_manual(values=colour_vec) + ggtitle(i)
}

ggRef + ggNew[["SDY296"]] + ggNew[["SDY301"]]
ggsave("../results/ns_umap_GMM.png",width=8,height=2.2)


# ---------------------------------
# Save cluster assignments
# ---------------------------------
clusts<-bind_rows(embed,embedNew[["SDY296"]],embedNew[["SDY301"]])
write.csv(clusts,"../results/GMM_k_3_D_2_clusters.csv")


# ---------------------------------
# Cluster characterisations
# ---------------------------------

# Could add Heatmap profiles of ns data for SDY96 and SDY301

# --------------------------------------------------------------------------------------------------------------------
# hai_identify_protected.R
# 	- having identified clusters then asses which clusters associate with protection by antibody readout
# 	- finds which Nanostring + demographic data defined clusters ultimately had antibody reponse post vaccine

options(tibble.print_max = 100)
options(tibble.print_width = Inf)

hai<-list()

hai[["SDY180"]] <-
read_csv("../data/sdy180/sdy180-dr47_tab/hai_result.csv") 

hai[["SDY296"]]<-
read_csv("../data/sdy296/resultfiles/hai_result.csv") 

hai[["SDY301"]]<-
  read_csv("../data/sdy301/resultfiles/hai_result.csv") 

# ---------------------------------
# Reduce to required fields
# ---------------------------------

fields<-c(grep("ACCESSION",names(hai[["SDY180"]]),value=TRUE),"VALUE_PREFERRED","STUDY_TIME_COLLECTED")

hai_df <- bind_rows(lapply(hai, function(df) df[, fields]))%>%select(-REPOSITORY_ACCESSION)

# ---------------------------------
# Read in clusters from ns UMAP GMM
# ---------------------------------
clusts<-read.csv("../results/GMM_k_3_D_2_clusters.csv",row.names=1)

# ---------------------------------
# Map sample id to subject
# ---------------------------------
clusts<-clusts%>%rename(SUBJECT_ACCESSION=Subject.Accession)
length(unique(hai_df$EXPSAMPLE_ACCESSION))

clust_hai<-inner_join(hai_df,clusts[,c("SUBJECT_ACCESSION","Group","Gender","Subject.Age")],by="SUBJECT_ACCESSION")

# ---------------------------------
# Plot cluster vs antibody response
# ---------------------------------
colour_vec<-c("#32CD32","#F37FB8","#409388","#CE8BAE","#B23648",
                  "#ADD8E6","#D46E7E","#7E486B","#79AA7A","#FFEC2B",
                  "#8D5B96","#E41A1C","#00B4F0","#3A85A8","#488846",
                  "#BD6253","#46A169","#EB7AA9","#C4625D","#D8B62E",
                  "#d6c624","#77777C","#4F6FA1","#E1712E","#A65628",
                  "#B392A3","#E984B9","#F2E631","#999999")
names(colour_vec) <- paste0("C",seq(from=1,to=length(colour_vec),by=1))

clust_hai$Group<-paste0("C",clust_hai$Group)
clust_hai$grp_time<-paste0(clust_hai$Group,clust_hai)

ggplot(clust_hai,aes(x=Group,y=VALUE_PREFERRED, fill=Group)) + geom_jitter(width = 0.1,size=0.2) + 
facet_wrap(~STUDY_TIME_COLLECTED) + geom_violin(width = 0.5,alpha = 0.7) +
ylab("HAI") + xlab("Cluster defined by UMAP-GMM on Nanostring data")  + scale_fill_manual(values=colour_vec) +
labs(fill="Cluster") + theme(legend.position="none") + stat_compare_means()
ggsave("../results/cluster_HAI.png",width=5,height=5)

# Day 28 only
clust_hai%>%filter(STUDY_TIME_COLLECTED==28)%>%
ggplot(.,aes(x=Group,y=VALUE_PREFERRED,fill=Group)) + geom_jitter(width = 0.1,size=0.2) + 
geom_violin(width = 0.5,alpha = 0.7) +
ylab("HAI") + xlab("Cluster defined by UMAP-GMM on Nanostring data") + scale_fill_manual(values=colour_vec) +
labs(fill="Cluster") + theme(legend.position="none") + stat_compare_means()
ggsave("../results/cluster_HAI_Day28.png",width=4,height=3)

# --------------------------------------------------------------------------------------------------------------------
# Read in new data to demonstrate prediction on single new case
# ---------------------------------
# Read in new data
# ---------------------------------
modelDir<-"../models/"
cluster_mod<-list.files(modelDir,pattern="GMM*")
umap_mod<-list.files(modelDir,pattern="UMAP*")

# ---------------------------------
# Read in new data
# ---------------------------------

new_sample<-read.csv("../data/randomNew.csv",row.names=1)

# ---------------------------------
# Load the transform and clustering model 
# ---------------------------------
umap_model<-uwot::load_uwot(file = paste0(modelDir,umap_mod))
gmm_model<-readRDS(paste0(modelDir,cluster_mod))

# ---------------------------------
# Apply models to get cluster assignment
# ---------------------------------
set.seed(10)
# Had to replicate row due to requirements for uwot::umap_transform
new_sample<-rbind(new_sample,new_sample)

UMAP_new <- uwot::umap_transform(new_sample,umap_model)
pred_new <- predict.Mclust(newdata = UMAP_new,object = gmm_model) # Fit model 

names(pred_new$classification)<-rownames(UMAP_new)

embed_new <- data.frame(X1=UMAP_new[1,1],X2=UMAP_new[1,2],Group=as.character(pred_new$classification[1]))
# ---------------------------------
# Plot new on initial embedding
# ---------------------------------
# Colour Palette ---------------------------------------------------
pal<-list()
pal$Group<-c("#32CD32","#F37FB8","#409388","#CE8BAE","#B23648",
                  "#ADD8E6","#D46E7E","#7E486B","#79AA7A","#FFEC2B",
                  "#8D5B96","#E41A1C","#00B4F0","#3A85A8","#488846",
                  "#BD6253","#46A169","#EB7AA9","#C4625D","#D8B62E",
                  "#d6c624","#77777C","#4F6FA1","#E1712E","#A65628",
                  "#B392A3","#E984B9","#F2E631","#999999")
pal$Study<-c(SDY180="brown",SDY296="purple", SDY301="grey")
pal$Arm<-c(ARM773="#E984B9",ARM776="#F2E631",ARM779="#999999")
pal$Time<-c("Day_-7"="#C4625D","Day_0"="#D8B62E")
pal$Gender<-c(Female="red",Male="blue")
names(pal$Group) <- seq(from=1,to=length(pal$Group),by=1)
# ------------------------------------------------------------------

# Add original meta to enable cluster labelling
expMeta<-"../data/sdy180/resultfiles/sdy180-dr47_subject_2_gene_expression_result.txt"
ns.meta<-read.delim(expMeta)

embed<-data.frame(umap_model$embedding,stringsAsFactors = F)

embed<-embed%>%rownames_to_column("Expsample.Accession")%>%left_join(.,ns.meta,by="Expsample.Accession")%>%
  mutate(Group=factor(gmm_model$classification))%>%mutate(Study.Time.Collected=paste0("Day_",Study.Time.Collected))

ggRef<-ggplot(embed,aes_string(x="X1",y="X2", colour="Group"))+
      geom_point(size=3,alpha=0.7)+
      theme_bw() +
      labs(color="GMM cluster")+
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(), legend.position="none")+
      xlab("UMAP1") + ylab("UMAP2")+
      scale_color_manual(values=pal$Group) + ggtitle("Reference UMAP+Clustering")

ggRefNew<-ggplot(embed,aes_string(x="X1",y="X2", colour="Group"))+
      geom_point(size=3,alpha=0.7)+
      theme_bw() +
      labs(color="GMM cluster")+
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank())+
      xlab("UMAP1") + ylab("UMAP2")+
      scale_color_manual(values=pal$Group) + ggtitle(paste0("New sample predicted cluster = ",pred_new$classification[1])) +
      geom_point(data=embed_new,size=7,alpha=0.7)+
      annotate("text", label = rownames(embed_new),x=embed_new$X1,y=embed_new$X2)


ggRef + ggRefNew
ggsave("../results/new_sample_umap_GMM.png",width=7,height=2.5)

# --------------------------------------------------------------------------------------------------------------------
# byCluster_profiles.R <MAYBE DONT RUN?>
# 	- taking assigned cluster labels and profiling these in other datasets having the same individuals
# 	- this can identify other readouts that are correlates of the protective cluster
# 	- in this case we look at RNA expression

clusts<-read.csv("../results/GMM_k_3_D_2_clusters.csv",row.names=1)

study<-c("SDY296","SDY301")
dataFiles <- c("../data/sdy296/resultfiles/rna_sequencing_result/SDY296_EXP13760_RNA_seq.703270.tsv",
  "../data/sdy301/resultfiles/rna_sequencing_result/SDY301_EXP13728_RNA_seq.703279.tsv")

metaFiles<-c(
  "../data/sdy296/resultfiles/sdy296-dr47_subject_2_rna_sequencing_result.txt",
  "../data/sdy301/resultfiles/sdy301-dr47_subject_2_rna_sequencing_result.txt")
# ---------------------------------
# Read Data
# ---------------------------------
gex.d<-lapply(dataFiles,read.delim)
names(gex.d)<-study

names(gex.d[[1]])[1]<-"ENSEMBL_ID"
gex<-inner_join(gex.d[[1]],gex.d[[2]], by="ENSEMBL_ID")%>%select(c("ENSEMBL_ID","GENE_SYMBOL",grep("^ES",names(.),value=TRUE)))

meta.d<-lapply(metaFiles,read.delim)
meta.d[[2]]$Planned.Visit.Name<-as.character(meta.d[[2]]$Planned.Visit.Name)
meta.d<-bind_rows(meta.d[[1]],meta.d[[2]])
meta.d<-meta.d%>%select(.,c("Subject.Accession","Gender","Subject.Age","ARM.Accession","Study.Accession","Expsample.Accession","Study.Time.Collected"))%>%
rename(id=Expsample.Accession)
# ---------------------------------
# Join meta
# ---------------------------------
clust_gex<-inner_join(meta.d,clusts[,c("Subject.Accession","Group","Gender","Subject.Age")],by="Subject.Accession")

# remove sample not in metadata
s<-grep("^ES",names(gex)[!names(gex)%in%meta.d$id],value=TRUE)
gex<-gex%>%select(-s)

# ---------------------------------
# Annotation data
# ---------------------------------
clust_gex<-clust_gex[clust_gex$id%in%names(gex),]
clust_gex<-distinct(clust_gex)

# remove if multiple samples from same subject 
s<-names(table(clust_gex$id))[table(clust_gex$id)>1] 
gex<-gex%>%select(-s)

clust_gex<-clust_gex[clust_gex$id%in%names(gex),]
keep<-names(gex)%in%clust_gex$id
keep[1:2]<-TRUE
gex<-gex[,keep]

# order annotation and gex by cluster
clust_gex<-clust_gex[order(clust_gex$Group),]
gex<-gex[,c("ENSEMBL_ID","GENE_SYMBOL",clust_gex$id)]
clust_gex$Group<-paste0("C",clust_gex$Group)
clust_gex$Study.Time.Collected<-paste0("Day_",clust_gex$Study.Time.Collected)
# ---------------------------------
# Gene selection
# ---------------------------------
# Could be marker genes for cluster identified as "protective"
n=100 # select most variable genes
d<-gex[,clust_gex$id]
rownames(d)<-gex$ENSEMBL_ID

gene_variance <- apply(d, 1, var)

# Sort the genes based on their variability
gene_variance <- gene_variance[order(gene_variance,decreasing=TRUE)]
#
sel<-names(gene_variance)[1:n]

# ---------------------------------
# Heatmap by cluster
# ---------------------------------
# Colour Palette ---------------------------------------------------
pal<-list()
pal$Group<-c("#32CD32","#F37FB8","#409388","#CE8BAE","#B23648",
                  "#ADD8E6","#D46E7E","#7E486B","#79AA7A","#FFEC2B",
                  "#8D5B96","#E41A1C","#00B4F0","#3A85A8","#488846",
                  "#BD6253","#46A169","#EB7AA9","#C4625D","#D8B62E",
                  "#d6c624","#77777C","#4F6FA1","#E1712E","#A65628",
                  "#B392A3","#E984B9","#F2E631","#999999")
pal$Study<-c(SDY180="brown",SDY296="purple", SDY301="grey")
pal$Arm<-c(ARM773="#E984B9",ARM776="#F2E631",ARM779="#999999",ARM2102="#8D5B96",ARM2107="#77777C")
pal$Time<-c("Day_-7"="#C4625D","Day_0"="#D8B62E","Day_1"="#00B4F0","Day_7"="#EB7AA9")
pal$Gender<-c(Female="red",Male="blue")
names(pal$Group) <- paste0("C",seq(from=1,to=length(pal$Group),by=1))
# ------------------------------------------------------------------

# need to set annotation colour palette
annotation_colors <- list(#Subject.Accession=palette.top$highLevelTypePal,
                          ARM.Accession=pal$Arm,
                          Study.Accession=pal$Study,
                          Study.Time.Collected=pal$Time,
                          Group=pal$Group)

row_ha = rowAnnotation(count = anno_boxplot(my_data,which="row",size = unit(0.5, "mm"), width = unit(1, "cm"),box_width = 0.3))

ha = HeatmapAnnotation(df = clust_gex[,c("Subject.Accession","ARM.Accession","Study.Accession","Study.Time.Collected","Group")],
                       show_annotation_name = TRUE,
                       col = annotation_colors,
                       simple_anno_size = unit(0.3, "cm"),
                       show_legend = c(FALSE,TRUE,TRUE,TRUE,TRUE))

# Expression data
my_data <- as.matrix(gex[,clust_gex$id])
rownames(my_data)<-gex$ENSEMBL_ID
my_data <- my_data[intersect(sel,rownames(my_data)),]
my_data<-t(scale(t(my_data)))
# Heatmap
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
# pdf(file = "../results/mvg_rna_heat.pdf",
#     width = 6,
#     height = 10, useDingbats = F)
# Heatmap(
#   my_data,
#   col = col_fun,
#   cluster_rows = TRUE,
#   cluster_columns = FALSE,
#   column_order = NULL,
#   show_row_dend = FALSE,
#   show_column_dend = FALSE,
#   show_row_names = TRUE,
#   show_column_names = FALSE,
#   bottom_annotation = NULL,
#   row_names_gp = gpar(fontsize = 6),
#   heatmap_legend_param = list(title = ""),
#   top_annotation = ha,
#   left_annotation = row_ha
# )
# dev.off()



	



