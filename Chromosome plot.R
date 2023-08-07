library(Gviz)
library(GenomicRanges)
data(geneModels)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library("org.Mm.eg.db")
library(rtracklayer)
geneRanges = import.gff("gencode.vM25.annotation.gtf")

itrack <- IdeogramTrack(genome = "mm10", chromosome = "chr19")

plotTracks(list(itrack),from = 10950000, to = 11650000)



geneRanges = geneRanges[grep("Ms4a",geneRanges$gene_name)]
geneRanges = geneRanges[geneRanges$type=="gene"]
geneRanges <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
geneRanges$symbol <- mapIds(org.Mm.eg.db,
                            keys = geneRanges$gene_id,
                            column = "SYMBOL",
                            keytype = "ENTREZID",
                            multiVals = "first")
geneRanges$symbol = geneRanges$gene_name


grtrack <- GeneRegionTrack(geneRanges, genome = "mm10",
                           chromosome = "chr19", name = "gene_name", 
                           transcriptAnnotation = "symbol",
                           background.panel = "white",
                           background.title = "darkblue")
pdf("Ch19_Ms4a_cluster.pdf", width = 10, height = 3)
plotTracks(list(itrack,grtrack),from = 10930000, to = 11650000)
dev.off()
