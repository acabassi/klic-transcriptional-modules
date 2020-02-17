################################ Transcriptional module discovery ##################################
######################################### Harbison data ############################################

# We wish to cluster genes, to identify those that are functionally related.
#
# We consider the example given in Section 3.3 of *Kirk, P., Griffin, J. E., Savage, R. S.,
# Ghahramani, Z., & Wild, D. L. (2012). Bayesian correlated clustering to integrate multiple
# datasets. Bioinformatics (Oxford, England), 28(24), 3290–3297.*

rm(list = ls())

set.seed(1)

library(BHC)
library(cluster)
library(ComplexHeatmap)
library(klic)
library(pheatmap)

########################################## Load the data ###########################################

harbisonDataNumeric <- read.table(
  "data/harbison_ppi_marina.csv",
  sep = ",",
  header = T,
  row.names = 1
)

harbisonData    <- harbisonDataNumeric
harbisonData[]  <- lapply(harbisonDataNumeric, factor)

########################################## Save data subsets #######################################

harbisonDataNumeric <-
  harbisonDataNumeric[, which(colSums(harbisonDataNumeric) > 0)]

# Remove string "-A" which causes problems when loading the csv files in matlab
odd_names <- rownames(harbisonDataNumeric)[grep("-A", rownames(harbisonDataNumeric))]
remove_part = function(string){
  new_string <- substr(string,1,nchar(string)-2)
  return(new_string)
}
replacement_names <- unlist(lapply(odd_names, remove_part))
to_be_replaced <- grep("-A", rownames(harbisonDataNumeric))
rownames(harbisonDataNumeric)[to_be_replaced] <- replacement_names

n_reps <- 200
pItem <- 0.95
N <- dim(harbisonDataNumeric)[1]

# for (i in 1:n_reps) {
#   items <- sample(N, floor(N * pItem))
#   write.csv(
#     harbisonDataNumeric[items,],
#     file = paste(
#       "GBNP/data/HarbisonDataSubsample",
#       i,
#       ".csv",
#       sep = ""
#     )
#   )
# }

######################### Do clustering in Matlab then load clusters into R ########################

# Co-clustering matrix
coClusteringMatrix <- matrix(0, N, N)
rownames(coClusteringMatrix) <-
  colnames(coClusteringMatrix) <- rownames(harbisonDataNumeric)
normalisingMatrix <- matrix(0, N, N)
rownames(normalisingMatrix) <-
  colnames(normalisingMatrix) <- rownames(harbisonDataNumeric)

for (i in 1:n_reps) {
  clusters <- read.csv(
    paste(
      "GBNP/data/HarbisonDataSubsample",
      i,
      "_clusters.csv",
      sep = ""
    ), row.names = 1,
  )
  inds <- match(rownames(clusters), rownames(coClusteringMatrix))
  inds <- inds[!is.na(inds)]
  maxK <- length(table(clusters))
  for (k in 1:maxK) {
    # Update counts of number of times each pair has been put in the same cluster
    coClusteringMatrix[inds, inds] <-
      coClusteringMatrix[inds, inds] +
      crossprod(t(as.numeric(clusters[rownames(clusters) %in% rownames(coClusteringMatrix),] == k)))
  }
  normalisingMatrix <-
    normalisingMatrix + crossprod(t(c(1:N) %in% inds))
}

consensusMatrix <- coClusteringMatrix / normalisingMatrix

pheatmap(consensusMatrix)

### Global clusters ###

clusters <- read.csv(
  paste(
    "GBNP/data/harbison_ppi_marina_clusters.csv",
    sep = ""
  ), header = FALSE, row.names = 1,
)

harbClusters <- clusters[,1]
rownames(clusters)[which(rownames(clusters) != rownames(harbisonDataNumeric))]
rownames(harbisonDataNumeric)[which(rownames(harbisonDataNumeric) != rownames(clusters))]
names(harbClusters) <- rownames(clusters)
rownames(harbisonDataNumeric)[to_be_replaced] <- odd_names

save(harbClusters, consensusMatrix, file = "data/harbisonClusters_GBNP.RData")
table(harbClusters)
sortedClusters <- sort(harbClusters, index.return = T)

col_bw = c("white", "black")
myBlues <-
  colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)

Hdata <-
  Heatmap(
    as.matrix(harbisonDataNumeric)[sortedClusters$ix, ],
    cluster_rows = F,
    cluster_columns = T,
    row_split = sortedClusters$x,
    show_row_names = FALSE,
    show_column_names = FALSE,
    name = "Data",
    col = col_bw,
    width = unit(5, "cm"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 28),
      labels_gp = gpar(fontsize = 28),
      direction = "horizontal",
      nrow = 1
    )
  )

Hclusters <-
  Heatmap(
    as.matrix(as.factor(sortedClusters$x)),
    show_row_names = FALSE,
    name = "Clusters",
    width = unit(.4, "cm"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 28),
      labels_gp = gpar(fontsize = 12),
      direction = "horizontal",
      nrow = 2
    )
  )

HcoClust <-
  Heatmap(
    consensusMatrix[sortedClusters$ix, sortedClusters$ix],
    cluster_rows = F,
    cluster_columns = F,
    row_split = sortedClusters$x,
    show_row_names = FALSE,
    show_column_names = F,
    col = myBlues,
    name = "Similarities",
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 28),
      labels_gp = gpar(fontsize = 28),
      legend_direction = "horizontal"
    )
  )
png("figures/harbison_initialClusters_GBNP_consMatrix.png",
    height = 730,
    width = 800
)
draw(
  Hdata + Hclusters + HcoClust,
  merge_legend = TRUE,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()

ccc <- copheneticCorrelation(consensusMatrix)
save(ccc, file = "cophenetic-correlation-coefficients/harbisonData_GBNP.RData")
