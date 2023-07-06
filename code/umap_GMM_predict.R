library(tidyverse)
library(uwot)
library(patchwork)
library(mclust)

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
