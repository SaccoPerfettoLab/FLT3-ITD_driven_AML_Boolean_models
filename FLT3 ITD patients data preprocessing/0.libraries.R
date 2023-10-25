
library(tidyverse)
library(maftools)


library(Rsubread)
library("DESeq2")

library(ggpubr)

library("ggdendro")

library(pheatmap)

library(Biostrings)

library(org.Hs.eg.db)

library(annotate)

# Custom function for plotPCA 

plotPCA1 <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  # object <- rld
  # intgroup = "condition"
  # ntop = 500
  # returnData = FALSE
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1){factor(apply(intgroup.df, 1, paste, collapse = ":"))}else {colData(object)[[intgroup]]}
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  library(ggrepel)
  ggplot(data = d, aes(x = PC1, y = PC2)) + 
    geom_point(size = 3, aes(color = group)) + 
    xlab(paste0("PC1: ", round(percentVar[1] *  100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    geom_text_repel(aes(label = name))+
    coord_fixed() 
}
