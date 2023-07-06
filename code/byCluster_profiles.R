# gene expression by cluster

library(tidyverse)
library(ComplexHeatmap)
library(patchwork)
library(RColorBrewer)
library(circlize)


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
pdf(file = "../results/example_rna_heat.pdf",
    width = 6,
    height = 10, useDingbats = F)
Heatmap(
  my_data,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = NULL,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  bottom_annotation = NULL,
  row_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(title = ""),
  top_annotation = ha,
  left_annotation = row_ha
)
dev.off()

