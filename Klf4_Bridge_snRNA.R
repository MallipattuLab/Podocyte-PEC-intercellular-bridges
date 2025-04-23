library(Seurat)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(harmony)
library(ggplot2)

setwd("Klf4_Datasets")

Klf4KO_151 <- Read10X_h5("Klf4/151_filtered_feature_bc_matrix.h5")
Klf4KO_151_Seurat_Object <- CreateSeuratObject(
  counts = Klf4KO_151,
  assay = "RNA",
  project = "Klf4_WT"
)

Klf4KO_161 <- Read10X_h5("Klf4/161_filtered_feature_bc_matrix.h5")
Klf4KO_161_Seurat_Object <- CreateSeuratObject(
  counts = Klf4KO_161,
  assay = "RNA",
  project = "Klf4_KO"
)

Klf4KO_33 <- Read10X_h5("Klf4/3_3_filtered_feature_bc_matrix.h5")
Klf4KO_33_Seurat_Object <- CreateSeuratObject(
  counts = Klf4KO_33,
  assay = "RNA",
  project = "Klf4_WT"
)

Klf4KO_34 <- Read10X_h5("Klf4/3_4_filtered_feature_bc_matrix.h5")
Klf4KO_34_Seurat_Object <- CreateSeuratObject(
  counts = Klf4KO_34,
  assay = "RNA",
  project = "Klf4_WT"
)


Klf4KO_41 <- Read10X_h5("Klf4/4_1_filtered_feature_bc_matrix.h5")
Klf4KO_41_Seurat_Object <- CreateSeuratObject(
  counts = Klf4KO_41,
  assay = "RNA",
  project = "Klf4_KO"
)

Klf4KO_42 <- Read10X_h5("Klf4/4_2_filtered_feature_bc_matrix.h5")
Klf4KO_42_Seurat_Object <- CreateSeuratObject(
  counts = Klf4KO_42,
  assay = "RNA",
  project = "Klf4_KO"
)

Klf4KO_151_Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Klf4KO_151_Seurat_Object, pattern = "^mt-")
Klf4KO_161_Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Klf4KO_161_Seurat_Object, pattern = "^mt-")
Klf4KO_33_Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Klf4KO_33_Seurat_Object, pattern = "^mt-")
Klf4KO_34_Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Klf4KO_34_Seurat_Object, pattern = "^mt-")
Klf4KO_41_Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Klf4KO_41_Seurat_Object, pattern = "^mt-")
Klf4KO_42_Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Klf4KO_42_Seurat_Object, pattern = "^mt-")

Klf4_Datasets.list <- lapply(
  X = c(
    Klf4KO_151_Seurat_Object,
    Klf4KO_161_Seurat_Object,
    Klf4KO_33_Seurat_Object,
    Klf4KO_34_Seurat_Object,
    Klf4KO_41_Seurat_Object,
    Klf4KO_42_Seurat_Object
    ), FUN = function(x) {
        
          x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)
          
          x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
        })


features <- SelectIntegrationFeatures(object.list = Klf4_Datasets.list)
Klf4_Datasets.anchors <- FindIntegrationAnchors(object.list = Klf4_Datasets.list, anchor.features = features)
Klf4_Datasets_combined <- IntegrateData(anchorset = Klf4_Datasets.anchors)

saveRDS(Klf4_Datasets_combined, file = "Klf4_Datasets_combined.rds")
Klf4_Datasets_combined <- readRDS("Klf4_Datasets_combined.rds")

Klf4_Datasets_combined$Sample <- Klf4_Datasets_combined$orig.ident

all.genes <- rownames(Klf4_Datasets_combined)
Klf4_Datasets_combined <- ScaleData(Klf4_Datasets_combined, features = all.genes)

Klf4_Datasets_combined <- RunPCA(Klf4_Datasets_combined, features = VariableFeatures(object = Klf4_Datasets_combined))

print(Klf4_Datasets_combined[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Klf4_Datasets_combined, dims = 1:2, reduction = "pca")

ElbowPlot(Klf4_Datasets_combined, ndims = 40)


Klf4_Datasets_combined <- RunHarmony(Klf4_Datasets_combined, "Sample", project.dim = F)
ElbowPlot(Klf4_Datasets_combined, ndims = 40, reduction = "harmony")

Klf4_Datasets_combined <- FindNeighbors(Klf4_Datasets_combined, dims = 1:23, reduction = "harmony")
Klf4_Datasets_combined <- FindClusters(Klf4_Datasets_combined, resolution = 0.8, reduction = "harmony")

Klf4_Datasets_combined <- RunUMAP(Klf4_Datasets_combined, reduction = "harmony", dims = 1:23, return.model = TRUE)

DimPlot(Klf4_Datasets_combined, reduction = "umap", label = TRUE, raster=FALSE)
DimPlot(Klf4_Datasets_combined, reduction = "umap", split.by = "orig.ident", label = TRUE, raster=FALSE) + NoLegend()

saveRDS(Klf4_Datasets_combined, file = "Klf4_Datasets_clustered.rds")
Klf4_Datasets_combined <- readRDS("Klf4_Datasets_clustered.rds")
# Klf4_Datasets_combined <- readRDS("TG26_clustered.rds")

Idents(Klf4_Datasets_combined) <- Klf4_Datasets_combined$ClusterName
Klf4_Datasets_combined.markers <- FindAllMarkers(Klf4_Datasets_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- Klf4_Datasets_combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

DotPlot(Klf4_Datasets_combined, features = c(
  "Sox2ot", 
  "Lrp2",
  "Slc25a36",
  "Slc27a2",
  "Dock10",
  "Calcrl", 
  "Umod", 
  "Fth1", 
  "Slc22a30",
  "Trpm7", 
  "Cfh", 
  "Slc4a9",  
  "Mki67", 
  "Wt1",
  "Frmpd4", 
  "Cd74", 
  "Wdr17",
  "Dmxl1", 
  "Meis2",
  "Akap12",
  "Arl15", 
  "Clu",  
  "Ebf1",
  "Thsd4"
),
        assay = "RNA", 
)+RotatedAxis()


new.cluster.ids <- c("PTS3-1",    #0
                     "PTS3-2",     #1
                     "PTS1-1",   #2 
                     "PTS1-2",    #3
                     "PTS1_S2",    #4
                     "CTAL",       #5
                     "PEC",    #6 
                     "EC1-1",       #7
                     "Tcell",     #8
                     "DCT-CNT-1",    #9
                     "CTAL",       #10
                     "DCT-CNT-2",    #11
                     "Per_Fib",        #12
                     "DCT",  #13
                     "ICB",   #14
                     "CNT",        #15
                     "Pod",  #16
                     "Prolif",   #17
                     "ICA",        #18
                     "Macrophage",  #19
                     "EC1-2",   #20
                     "Per",   #21
                     "EC2",    #22
                     "Uro"   #23
)
names(new.cluster.ids) <- levels(Klf4_Datasets_combined)
Klf4_Datasets_combined <- RenameIdents(Klf4_Datasets_combined, new.cluster.ids)


saveRDS(Klf4_Datasets_combined, file = "Klf4_Datasets_labeled.rds")
Klf4_Datasets_combined <- readRDS("Klf4_Datasets_labeled_version0.rds")
Klf4_Datasets_combined$ClusterName <- Idents(Klf4_Datasets_combined)

Idents(Klf4_Datasets_combined) <- Klf4_Datasets_combined$ClusterName
DimPlot(Klf4_Datasets_combined, reduction = "umap", label = TRUE, raster=FALSE)
DimPlot(Klf4_Datasets_combined, reduction = "umap", split.by = "Sample", label = TRUE, raster=FALSE) + NoLegend()

Klf4_Datasets_combined$Group <- Klf4_Datasets_combined$Sample
Klf4_Datasets_combined$Group[which(Klf4_Datasets_combined$Group %in% c("Klf4_WT"))] <- "WT"
Klf4_Datasets_combined$Group[which(Klf4_Datasets_combined$Group %in% c("Klf4_KO"))] <- "Injury"

####################### Subcluster ######################################
Idents(Klf4_Datasets_combined ) <- Klf4_Datasets_combined$ClusterName

Klf4_Datasets_combined <- FindSubCluster(
  Klf4_Datasets_combined,
  c("PEC"),
  graph.name = 'integrated_nn',
  subcluster.name = "PEC",
  resolution = 0.4,
  algorithm = 1
)

Idents(Klf4_Datasets_combined) <- Klf4_Datasets_combined$PEC

Klf4_Datasets_combined.markers <- FindAllMarkers(Klf4_Datasets_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

results <- list()
clusters <- unique(Klf4_Datasets_combined.markers$cluster)
for(cluster in clusters) {
  markers <- Klf4_Datasets_combined.markers[which(Klf4_Datasets_combined.markers$cluster == cluster),]
  # markers$gene <- rownames(markers)
  cluster <- gsub("/", "_", cluster)
  results[[cluster]] <- markers
}

new.cluster.ids <- c(
  "PTS3-1",
  "PTS1-1",
  "PTS3-2",
  "Per_Fib",
  "CTAL",
  "ICA",
  "DCT",
  "EC2",
  "ICB",
  "PEC_1",
  "Mito-enriched",
  "PTS1-2",
  "EC1-2",
  "PC",
  "DCT-CNT",
  "EC1-1",
  "PEC_2",
  "PTS3-3",
  "Macrophage",
  "PEC_3",
  "Pod",
  "PTS3-PEC",
  "PTS1-PEC",
  "Prolif",
  "Uro",
  "PEC_0",
  "Quienscient_PEC",
  "MD",
  "Per",
  "PTS1-1",
  "PTS3-PEC"  
)
names(new.cluster.ids) <- levels(Klf4_Datasets_combined)
Klf4_Datasets_combined <- RenameIdents(Klf4_Datasets_combined, new.cluster.ids)

Klf4_Datasets_combined$PEC_labeled <- Idents(Klf4_Datasets_combined)

Idents(Klf4_Datasets_combined) <- Klf4_Datasets_combined$PEC_labeled

cluster_composition <- data.frame(cluster = Klf4_Datasets_combined$PEC_labeled,
                                  sample = Klf4_Datasets_combined$Group
)

data <- as.data.frame(table(cluster_composition))
data <- data[order(data$cluster, decreasing = TRUE),]

g <- ggplot(data, aes(fill=sample, y=Freq, x=cluster))+ geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle=45, hjust=1))
g

for(Sample in sort(unique(Klf4_Datasets_combined$Sample))){
  data[which(data$sample %in% Sample),]$Freq <- data[which(data$sample %in% Sample),]$Freq/sum(data[which(data$sample %in% Sample),]$Freq)
}

g2 <- ggplot(data, aes(fill=sample, y=Freq, x=cluster)) +
  geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(angle=45, hjust=1))
g2

group_composition <- data.frame(group = Klf4_Datasets_combined$Group,
                                  sample = Klf4_Datasets_combined$Sample
)

data <- as.data.frame(table(group_composition))
data <- data[order(data$sample, decreasing = TRUE),]
data <- data[which(data$Freq != 0),]

g <- ggplot(data, aes(fill=group, y=Freq, x=sample))+ geom_bar(position="stack", stat="identity") + 
  geom_text(data=data,aes(x=sample,y=Freq,label=Freq),vjust=0) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
g

########################## Pathway Analysis ################################################
Idents(Klf4_Subset) <- Klf4_Subset$ClusterName
VlnPlot(Klf4_Subset, features = "Birc3", idents = c("PEC", "Pod"), assay = "RNA", split.by = "Sample",split.plot = TRUE, adjust = 8)

Klf4_Subset.markers <- FindAllMarkers(Klf4_Subset, assay = "SCT", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

results <- list()
clusters <- unique(Klf4_Subset.markers$cluster)
for(cluster in clusters) {
  markers <- Klf4_Subset.markers[which(Klf4_Subset.markers$cluster == cluster),]
  # markers$gene <- rownames(markers)
  cluster <- gsub("/", "_", cluster)
  results[[cluster]] <- markers
}
openxlsx::write.xlsx(results, "PEC_Subset_DEGs_Klf4_Only_SCT_UPDOWN.xlsx")

Klf4_Subset$Sample <- factor(x = Klf4_Subset$Sample, levels = c("Klf4_WT", "Klf4_KO"))

Klf4_Subset.markers <- FindAllMarkers(Klf4_Subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

results <- list()
clusters <- unique(Idents(Klf4_Subset))
clusters <- c("PEC", "Pod")
for(cluster in clusters) {
  markers <- FindMarkers(Klf4_Subset, ident.1 = cluster, min.pct = 0.25, logfc.threshold = 0.25)
  markers$gene <- rownames(markers)
  cluster <- gsub("/", "_", cluster)
  results[[cluster]] <- markers
}



library(ChIPseeker)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)
OrgDb = org.Mm.eg.db
library(ggplot2)
require(DOSE)
# We will lose some genes here because not all IDs will be converted
df <-  results$PEC[results$PEC$avg_log2FC > 0,]
ids<-bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=OrgDb) # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[(df$gene %in% dedup_ids$SYMBOL),]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
gene_list <- df2$avg_log2FC

# Name vector with ENTREZ ids
names(gene_list) <- df2$Y

# omit any NA values 
gene_list<-na.omit(gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results from df2
kegg_sig_genes_df = subset(df2, p_val_adj < 0.05)

# From significant results, we want to filter on log2fold change
genes <- kegg_sig_genes_df$avg_log2FC

# Name the vector with the CONVERTED ID!
names(genes) <- kegg_sig_genes_df$Y

# omit NA values
genes <- na.omit(genes)

kegg_organism = "mmu"
kk <- enrichKEGG(gene=names(gene_list), organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
kk_readable <- setReadable(kk, 'org.Mm.eg.db', 'ENTREZID')

barplot(kk, 
        showCategory = 15, 
        title = "KEGG Enriched Pathways SNP",
        font.size = 8)

dotplot(kk, 
        showCategory = 15, 
        title = "KEGG Enriched Pathways SNP",
        font.size = 8)
ggsave(filename = "Klf4_KEGG_dotplot_pureEnrichment.pdf",width = 8, height = 8)

# We will lose some genes here because not all IDs will be converted
df <-  results$Pod[results$Pod$avg_log2FC > 0,]
ids<-bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=OrgDb) # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[(df$gene %in% dedup_ids$SYMBOL),]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
gene_list <- df2$avg_log2FC

# Name vector with ENTREZ ids
names(gene_list) <- df2$Y

# omit any NA values 
gene_list<-na.omit(gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results from df2
kegg_sig_genes_df = subset(df2, p_val_adj < 0.05)

# From significant results, we want to filter on log2fold change
genes <- kegg_sig_genes_df$avg_log2FC

# Name the vector with the CONVERTED ID!
names(genes) <- kegg_sig_genes_df$Y

# omit NA values
genes <- na.omit(genes)

kegg_organism = "mmu"
kk <- enrichKEGG(gene=names(gene_list), organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
kk_readable <- setReadable(kk, 'org.Mm.eg.db', 'ENTREZID')

barplot(kk, 
        showCategory = 15, 
        title = "KEGG Enriched Pathways SNP",
        font.size = 8)

dotplot(kk, 
        showCategory = 15, 
        title = "KEGG Enriched Pathways SNP",
        font.size = 8)
ggsave(filename = "Klf4_KEGG_dotplot_pureEnrichment_Pod.pdf",width = 8, height = 8)


Idents(Klf4_Subset) <- Klf4_Subset$ClusterName
Klf4_Subset_PEC_Pod <- subset(Klf4_Subset, idents =c("PEC", "Pod") )
DefaultAssay(Klf4_Subset_PEC_Pod) <- "SCT"
# Klf4_Subset_PEC_Pod$ClusterSplited <- paste0(Klf4_Subset_PEC_Pod$ClusterName, "_", Klf4_Subset_PEC_Pod$Sample)

Klf4_Subset_PEC_Pod.bulk <- AverageExpression(Klf4_Subset_PEC_Pod, return.seurat = TRUE, group.by = c("ClusterSplited"))

Klf4_Subset_PEC_Pod.bulk <- ScaleData(Klf4_Subset_PEC_Pod.bulk, features = genes_in_pathway)

Idents(Klf4_Subset_PEC_Pod.bulk) <- factor(x = names(Idents(Klf4_Subset_PEC_Pod.bulk)), levels = c("PEC_Klf4_WT", "PEC_Klf4_KO","Pod_Klf4_WT","Pod_Klf4_KO"))
# Klf4_Subset_PEC_Pod.bulk$clusterSplited <- names(Idents(Klf4_Subset_PEC_Pod.bulk))
DoHeatmap(Klf4_Subset_PEC_Pod.bulk, features = unique(genes_in_pathway), assay = "SCT", 
          draw.lines = FALSE,
          slot = "scale.data")+
  scale_fill_gradient2(low="blue", high="red", 
                       midpoint=0, limits=range(c(-2,2.5)))


DoHeatmap(Klf4_Subset_PEC_Pod.bulk, features = c("Itgav","Birc3","Vcl","Itgb6","Itga1","Pdgfd","Col4a1","Mapk1","Pak2"), #Focal adhesion
          draw.lines = FALSE,
          assay = "SCT", slot = "scale.data")+
  scale_fill_gradient2(low="blue", high="red", 
                       midpoint=0, limits=range(c(-2,2.5)))
ggsave(filename = "Klf4_KEGG_Focal_adhesion.pdf",width = 5, height = 8)

DoHeatmap(Klf4_Subset_PEC_Pod.bulk, features = c("Nfkb1","Birc3","Malt1","Chuk","Plau","Ikbkb"), #NF-kappa B signaling pathway
          draw.lines = FALSE,
          assay = "SCT", slot = "scale.data")+
  scale_fill_gradient2(low="blue", high="red", 
                       midpoint=0, limits=range(c(-2,2.5)))
ggsave(filename = "Klf4_KEGG_NFKB_signaling.pdf",width = 5, height = 8)

DoHeatmap(Klf4_Subset_PEC_Pod.bulk, features = c("Vcam1","Nfkb1","Cxcl1","Creb5","Birc3","Birc2","Ikbkb"), #TNF signaling pathway
          draw.lines = FALSE,
          assay = "SCT", slot = "scale.data")+
  scale_fill_gradient2(low="blue", high="red", 
                       midpoint=0, limits=range(c(-2,2.5)))
ggsave(filename = "Klf4_KEGG_TNF_signaling.pdf",width = 5, height = 8)

DoHeatmap(Klf4_Subset_PEC_Pod.bulk, features = c("Itgb6","Nfkb1","Cdk6","Creb5","Pdgfd","Pdgfb","Itgav","Col4a1","Itga1","Ikbkb","Pkn2"), #PI3K-Akt signaling pathway
          draw.lines = FALSE,
          assay = "SCT", slot = "scale.data")+
  scale_fill_gradient2(low="blue", high="red", 
                       midpoint=0, limits=range(c(-2,2.5)))
ggsave(filename = "Klf4_KEGG_PI3K-Akt_signaling.pdf",width = 5, height = 8)

DoHeatmap(Klf4_Subset_PEC_Pod.bulk, features = c("Ntng2","Sema3e","Ephb1","Fyn","Srgap1","Sema3g","Robo2","Unc5c","Nck2","Srgap2","Pak1","Ptch1","Ablim1","Ryk"), #Axon guidance
          draw.lines = FALSE,
          assay = "SCT", slot = "scale.data")+
  scale_fill_gradient2(low="blue", high="red", 
                       midpoint=0, limits=range(c(-2,2.5)))
ggsave(filename = "Klf4_KEGG_Axon_guidance.pdf",width = 5, height = 8)

DoHeatmap(Klf4_Subset_PEC_Pod.bulk, features = c("Itgav","Rras2","Vcl","Itgb6","Itga1","Pdgfd","Myh9","Mapk1","Pak2"), #Regulation of actin cytoskeleton
          draw.lines = FALSE,
          assay = "SCT", slot = "scale.data")+
  scale_fill_gradient2(low="blue", high="red", 
                       midpoint=0, limits=range(c(-2,2.5)))
ggsave(filename = "Klf4_KEGG_Regulation_actin_cytoskeleton.pdf",width = 5, height = 8)

DoHeatmap(Klf4_Subset_PEC_Pod.bulk, features = c("Itgav","Vcl","Itgb6","Itga1","Agrn","Col4a1","Myh9","Sntb2","Daam"), #Cytoskeleton in muscle cells
          draw.lines = FALSE,
          assay = "SCT", slot = "scale.data")+
  scale_fill_gradient2(low="blue", high="red", 
                       midpoint=0, limits=range(c(-2,2.5)))
ggsave(filename = "Klf4_KEGG_Cytoskeleton_in_muscle_cells.pdf",width = 5, height = 8)

##################################### GSEA ################################

df <- FindMarkers(Klf4_Subset,ident.1 =c("PEC") , only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

df <- df[which(df$p_val_adj < 0.05),]
df <- df[order(df$avg_log2FC, decreasing = FALSE), ] 
df$gene <- rownames(df)

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(df$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb = org.Mm.eg.db)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)

dedup_ids = ids[!duplicated(ids$ENTREZID),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$gene %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$ENTREZID = dedup_ids$ENTREZID[match(df2$gene, dedup_ids$SYMBOL)]

df2 <- df2[is.finite(df2$avg_log2FC),]

# Create a vector of the gene unuiverse
gene_list <- df2$avg_log2FC

# Name vector with ENTREZ ids
names(gene_list) <- df2$ENTREZID

# omit any NA values 
gene_list<-na.omit(gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gse <-gseKEGG(geneList     = gene_list,
              organism     = "mmu",
              minGSSize    = 3,
              maxGSSize    = 800,
              pvalueCutoff = 0.05,
              pAdjustMethod = "none",
              keyType       = "ncbi-geneid")
gse_readable <- setReadable(gse, 'org.Mm.eg.db', 'ENTREZID')
gse_UP <- gse
gse_UP@result <- gse_UP@result[which(gse_UP@result$enrichmentScore > 0),]

gse_readable_UP <- setReadable(gse_UP, 'org.Mm.eg.db', 'ENTREZID')

dotplot(gse_readable_UP, showCategory=8, split=".sign") + facet_grid(.~.sign) + ggtitle("PEC")
ggsave(filename = "Klf4_KEGG_GSEA_dotplot.pdf",width = 8, height = 5)
options(enrichplot.colours = c("#327eba", "white", "#e06663"))
heatplot(gse_readable_UP, foldChange=gene_list, showCategory=2)+ ggtitle("PEC")
ggsave(filename = "Klf4_KEGG_GSEA_heatmap.pdf",width = 5, height = 2)


gse_readable_selected <- gse_readable
heatplot(gse_readable_selected, foldChange=gene_list, showCategory=16)+ ggtitle("PEC")
ggsave(filename = "Klf4_KEGG_GSEA.pdf",width = 5, height = 2)

############################# Klf4 Pod Subset ###################################

Idents(Klf4_Subset) <- Klf4_Subset$ClusterName
Pod_Subset <- subset(Klf4_Subset, idents ="Pod" )

DotPlot(Pod_Subset, assay = "SCT", features = "Birc3", split.by = "Sample")
ggsave(filename = "Klf4_Bric3_in_Pod.pdf",width = 4, height = 5)
