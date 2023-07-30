library(Gviz)
library(GenomicRanges)
data(cpgIslands)
class(cpgIslands)
chr <- as.character(unique(seqnames(cpgIslands)))
gen <- genome(cpgIslands)
atrack <- AnnotationTrack(cpgIslands, name = "CpG")
plotTracks(atrack)
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, atrack))
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
plotTracks(list(itrack, gtrack, atrack))
data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome = gen,
                           chromosome = chr, name = "Gene Model")
plotTracks(list(itrack, gtrack, atrack, grtrack))
plotTracks(list(itrack, gtrack, atrack, grtrack),
           from = 26700000, to = 26750000)
plotTracks(list(itrack, gtrack, atrack, grtrack),
           extend.left = 0.5, extend.right = 1000000)
plotTracks(list(itrack, gtrack, atrack, grtrack), 
           extend.left = 0.5, extend.right = 1000000, col = NULL)
library(BSgenome.Hsapiens.UCSC.hg38)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
?available.genomes
