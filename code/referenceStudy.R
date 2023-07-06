library(tidyverse)
library(uwot)
library(mclust)
library(ComplexHeatmap)
library(patchwork)
library(RColorBrewer)
library(circlize)

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

uwot::save_uwot(umap_model,file = paste0(modelDir,"UMAP_",D_value,"D_model_export.rds"))

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
      geom_point(size=1,alpha=0.7)+
      labs(x="UMAP1", y="UMAP2") +
      theme_bw() +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())#+
      #scale_color_manual(values=c("#00BFC4","#F8766D"))
    gg

embed<-embed%>%rownames_to_column("Expsample.Accession")%>%left_join(.,ns.meta,by="Expsample.Accession")%>%
  mutate(Group=factor(gmm_model$classification))%>%mutate(Study.Time.Collected=paste0("Day_",Study.Time.Collected))


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
ggsave("../results/SDY180_ns_umap.png",width=8,height=2.2)


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



