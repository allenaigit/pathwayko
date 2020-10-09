#' @import org.Mm.eg.db
#' @import Rgraphviz
#' @importFrom KEGGgraph parseKGML mergeGraphs getName getSubtype getType getKEGGedgeData KEGGpathway2Graph
#' @importFrom Biobase exprs annotation pData pData<- sampleNames<- sampleNames
#' @importFrom annotate getEG
#' @importFrom oligo read.celfiles rma
#' @importFrom AnnotationDbi select
#' @importFrom limma contrasts.fit eBayes lmFit makeContrasts topTable
#' @importFrom methods show new Summary
#' @importFrom stats pairwise.t.test pairwise.wilcox.test shapiro.test model.matrix median na.omit setNames
#' @importFrom utils read.csv untar write.csv
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics points legend text abline hist
#' @importFrom ggplot2 theme_bw theme element_text element_blank ggplot aes geom_violin geom_boxplot stat_summary geom_point
#' @importFrom gridExtra grid.arrange
#' @importFrom graph nodes edges nodes<- edgeDataDefaults<- edgeData edgeData<- 
#' @importFrom magick image_read image_write
#' @importFrom pROC roc coords plot.roc ci.sp roc.test
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils getTxtProgressBar setTxtProgressBar txtProgressBar
#' @importFrom parallel detectCores mclapply mccollect mcparallel
#' @importFrom stringr str_split
#' @importFrom igraph graph.data.frame igraph.to.graphNEL
#' @importFrom changepoint cpt.meanvar

#' @importFrom SPIA plotP spia
#' @importFrom ROntoTools setNodeWeights alphaMLG pe pDis
#' @importFrom safe getCmatrix safe safe.toptable
#' @importFrom fgsea fgsea
#' @importFrom PADOG padog
#' @importFrom GSA GSA

NULL