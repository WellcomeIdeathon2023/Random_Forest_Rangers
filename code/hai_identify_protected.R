library(tidyverse)
library(patchwork)
library(ggpubr)

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
facet_wrap(~STUDY_TIME_COLLECTED) + geom_violin(width = 0.5,alpha = 0.7, outlier.shape = NA) +
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

# Day 28 only facet by study
clust_hai%>%filter(STUDY_TIME_COLLECTED==28)%>%
ggplot(.,aes(x=Group,y=VALUE_PREFERRED,fill=Group)) + geom_jitter(width = 0.1,size=0.2) + 
geom_violin(width = 0.5,alpha = 0.7) + facet_wrap(~STUDY_ACCESSION) +
ylab("HAI") + xlab("Cluster defined by UMAP-GMM on reference data (SDY180)") + scale_fill_manual(values=colour_vec) +
labs(fill="Cluster") + theme(legend.position="none") + stat_compare_means()
ggsave("../results/cluster_HAI_Day28_byStudy.png",width=8,height=3)



