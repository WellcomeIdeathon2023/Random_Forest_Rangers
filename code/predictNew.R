library(tidyverse)
library(uwot)
library(patchwork)
library(mclust)

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






