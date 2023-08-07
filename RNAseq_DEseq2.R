library(DESeq2)
library(ggplot2)
SampleName = list.files(path = "star/",pattern = "R.bw")
SampleName = gsub("_sorted_R.bw","",SampleName)
Group = gsub("[0-9]$","",SampleName)
sample.meta = cbind.data.frame(SampleName,Group)
sample.meta = apply(sample.meta, 2, as.factor)
mtx = read.table("genes_expression_expected_count.tsv", sep = "\t", header = T, row.names = 1)
mtx = mtx[,-1]
df2 <- data.frame(apply(mtx, 2, function(x) round(x,0)))


dds <- DESeqDataSetFromMatrix(countData = df2,
                              colData = sample.meta,
                              design= ~ Group)
saveRDS(dds,"dds.rds")
##filtering
#keep <- rowSums(counts(dds)) >= 20
#dds <- dds[keep,]
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)

saveRDS(vsd,"vsd.rds")

library(pheatmap)
library(RColorBrewer)
library(viridis)
select <- unique(grep("Ms4a",rownames(ntd), value = T))
p = pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=T, color = brewer.pal(9, "OrRd"))
ggsave(plot = p, filename = "plots/Heat_map_Ms4as.jpeg", width = 3, height = 3, units = "in")


p = pheatmap(assay(ntd)[select,grep("WT|CL",colnames(assay(ntd)))], cluster_rows=T, show_rownames=T,
             cluster_cols=T, color = brewer.pal(9, "OrRd"))
ggsave(plot = p, filename = "plots/Heat_map_Ms4as_1.pdf", width = 3, height = 4, units = "in")

p = pheatmap(assay(ntd)[select,grep("WT|C6",colnames(assay(ntd)))], cluster_rows=T, show_rownames=T,
             cluster_cols=T, color = brewer.pal(9, "OrRd"))
ggsave(plot = p, filename = "plots/Heat_map_Ms4as_2.pdf", width = 3, height = 4, units = "in")


## keep all gene
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)

select <- unique(grep("Ms4a",rownames(ntd), value = T))
p1 = pheatmap(assay(ntd)[select,grep("WT|CL",colnames(assay(ntd)))], cluster_rows=T, show_rownames=T,
             cluster_cols=T, color = brewer.pal(9, "OrRd"))
ggsave(plot = p1, filename = "plots/Heat_map_Ms4as_WT_and_cluster.jpeg", width = 3, height = 5, units = "in")

p2 = pheatmap(assay(ntd)[select,grep("WT|C6",colnames(assay(ntd)))], cluster_rows=T, show_rownames=T,
              cluster_cols=T, color = brewer.pal(9, "OrRd"))
ggsave(plot = p2, filename = "plots/Heat_map_Ms4as_WT_and_6CKO.jpeg", width = 3, height = 5, units = "in")

pheatmap(assay(ntd)[select,grep("WT|C6",colnames(assay(ntd)))], cluster_rows=T, show_rownames=T,
         cluster_cols=T, color = viridis(100, option = "magma", direction = -1))
