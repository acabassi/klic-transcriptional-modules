################################ Transcriptional module discovery ##################################
######################################### Harbison data ############################################

# We wish to cluster genes, to identify those that are functionally related.
# 
# We consider the example given in Section 3.3 of *Kirk, P., Griffin, J. E., Savage, R. S., 
# Ghahramani, Z., & Wild, D. L. (2012). Bayesian correlated clustering to integrate multiple 
# datasets. Bioinformatics (Oxford, England), 28(24), 3290–3297.*

rm(list=ls())

set.seed(1)

library(BHC)
library(ComplexHeatmap)
library(klic)
library(pheatmap)

source("functions/getClusters.R")

########################################## Load the data ###########################################

harbisonDataNumeric <- read.table("data/harbison_ppi_marina.csv", 
                                  sep = ",", header = T, row.names = 1)

harbisonData    <- harbisonDataNumeric
harbisonData[]  <- lapply(harbisonDataNumeric, factor)

################################## Bayesian Hierarchical Clustering ################################

cat("Harbison data: Bayesian Hierarchical Clustering \n")

harbisonDataNumeric <- harbisonDataNumeric[,which(colSums(harbisonDataNumeric) > 0)]

save(harbisonDataNumeric, file = "data/harbisonData.RData")

n_reps <- 200
pItem <- 0.95
N <- dim(harbisonDataNumeric)[1]

# Co-clustering matrix
coClusteringMatrix <- matrix(0, N, N)
rownames(coClusteringMatrix) <- colnames(coClusteringMatrix) <- rownames(harbisonDataNumeric)

# Normalising matrix
normalisingMatrix <- matrix(0, N, N)
rownames(normalisingMatrix) <- colnames(normalisingMatrix) <- rownames(harbisonDataNumeric)

# Let's store all clusters so that we can plot them one-by-one later on
harbClusters <- matrix(NA, n_reps, N)
colnames(harbClusters) <- rownames(harbisonDataNumeric)

# CAUTION: This takes a long time to run. You can load previously output below.
# for(i in 1:n_reps){
# 
#   cat(i, "\n")
#   set.seed(i)
#   items <- sample(N, floor(N*pItem))
#   normalisingMatrix <- normalisingMatrix + crossprod(t(c(1:N)%in%items))
# 
#   bhc_harb <- bhc(harbisonDataNumeric[items,], itemLabels = rownames(harbisonDataNumeric)[items],
#                   randomised = FALSE, verbose = T)
#   harbClusters[i,items] <- getClusters(bhc_harb)
# 
#   maxK <- length(table(harbClusters))
#   for(k in 1:maxK){
#     # Update counts of number of times each pair has been put in the same cluster
#     coClusteringMatrix[items,items] <- coClusteringMatrix[items,items] +
#     crossprod(t(as.numeric(harbClusters[i,items]==k)))
#   }
# }
# 
# consensusMatrix <- coClusteringMatrix/normalisingMatrix
# 
# pheatmap(consensusMatrix)

# save(harbClusters,  file = paste("harbClusters_all_",pItem,"percentData_BHC.RData", sep = ""))
# load(paste("harbClusters_all_",pItem,"percentData_BHC.RData", sep = ""))

load("data/harbisonClusters_BHC.RData")

#################################### Plot intermediate clusters ####################################

# harbisonDataNumeric <- as.matrix(apply(harbisonDataNumeric, 2, as.factor))
# col_fun = c("white", "black")
# 
# n_plots <- 10
# for(i in 1:n_plots){
#     items <- which(!is.na(harbClusters[i,]))
#     sortedClusters <- sort(harbClusters[i,], index.return = T)
#     n_clusters <- max(sortedClusters$x)
#     rowGaps <- c()
#     for(j in 1:n_clusters){
#         rowGaps <- c(rowGaps, match(j, sortedClusters$x))
#     }
#     rowGaps <- rowGaps-1
#     Hdata <- Heatmap(harbisonDataNumeric[sortedClusters$ix,], cluster_rows = F,
#                      cluster_columns = T, row_split = sortedClusters$x,
#                      show_row_names = FALSE, name = "Data", col = col_fun)
#     Hclusters <- Heatmap(as.matrix(as.factor(sortedClusters$x)), show_row_names = FALSE,
#                          name = "Clusters")
#     png(paste("figures/harbisonData_initialClusters", i, "-", pItem,"percentData_BHC.png",sep=""),
#         height = 700, width = 800)
#     Hdata + Hclusters
#     dev.off()
# }

#################################### Clustering all observations ###################################

# bhc_harb <- bhc(harbisonDataNumeric, itemLabels = rownames(harbisonDataNumeric),
#                 randomised = FALSE, verbose = T)
# harbClusters <- getClusters(bhc_harb)
# 
# save(harbClusters, consensusMatrix, theta, file = "data/harbisonClusters_BHC.RData")
# 
# write.table(
#     as.data.frame(harbClusters),
#     paste("goto-scores/harbison_initialClusters_BHC.csv", sep = ""),
#     col.names = FALSE,
#     quote = TRUE
# )

# load("data/harbisonClusters_BHC.RData")

######################################## Plot final clusters #######################################

inds <- match(rownames(harbisonDataNumeric), names(harbClusters))
harbClusters <- harbClusters[inds]

inds <- match(rownames(harbisonDataNumeric), rownames(consensusMatrix))
sum(1-(inds == 1:551))
  
# Let's just plot the clusters:
table(harbClusters)
sortedClusters <- sort(harbClusters, index.return = T)

# Palettes
col_bw = c("white", "black")

cbPalette <-
  c(
    "#999999",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  )

# Plot as a heatmap
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
      labels_gp = gpar(fontsize = 28),
      direction = "horizontal",
      nrow = 1
    ),
    col = structure(cbPalette[1:5], names=c("1", "2", "3", "4", "5"))
  )

myBlues <-
    colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)
HcoClust <-
    Heatmap(
        consensusMatrix[sortedClusters$ix, sortedClusters$ix],
        cluster_rows = F,
        cluster_columns = F,
        col = myBlues,
        row_split = sortedClusters$x,
        column_split = sortedClusters$x,
        show_row_names = FALSE,
        show_column_names = FALSE,
        name = "Similarities",
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 28),
          labels_gp = gpar(fontsize = 28),
          legend_direction = "horizontal"
        )
    )

png("figures/harbison_initialClusters_BHC_consMatrix.png",
    height = 700,
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
save(ccc, file = "cophenetic-correlation-coefficients/harbisonData_BHC.RData")
 