library(tidyverse)
library(patchwork)

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

print(tmp, width = Inf)

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
clust_hai$Group<-paste0("C",clust_hai$Group)
clust_hai$grp_time<-paste0(clust_hai$Group,clust_hai)
ggplot(clust_hai,aes(x=Group,y=VALUE_PREFERRED)) + geom_jitter(width = 0.1,size=0.2) + 
facet_wrap(~STUDY_TIME_COLLECTED) + geom_violin(width = 0.5,alpha = 0.5, outlier.shape = NA) +
ylab("HAI") + xlab("Nanostring data defined cluster")
ggsave("../results/cluster_HAI.png",width=5,height=4.5)
