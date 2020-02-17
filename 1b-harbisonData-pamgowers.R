################################ Transcriptional module discovery ##################################
######################################### Harbison data ############################################

# We wish to cluster genes, to identify those that are functionally related.
#
# We consider the example given in Section 3.3 of *Kirk, P., Griffin, J. E., Savage, R. S.,
# Ghahramani, Z., & Wild, D. L. (2012). Bayesian correlated clustering to integrate multiple
# datasets. Bioinformatics (Oxford, England), 28(24), 3290–3297.*

rm(list = ls())

set.seed(1)
library(cluster)
library(pheatmap)
library(ComplexHeatmap)

########################################## Load the data ###########################################

harbisonDataNumeric <- read.table(
  "data/harbison_ppi_marina.csv",
  sep = ",",
  header = T,
  row.names = 1
)

################################## Partitioning Around Medoids ################################

cat("Harbison data: PAM with Gower's distance \n")

harbisonDataNumeric <- harbisonDataNumeric[,which(colSums(harbisonDataNumeric) > 0)]

distHarb <- daisy(harbisonDataNumeric, metric = c("gower"))
clusteringList <- clusterStatsList <- list(length = 20)

for(i in 2:20){
  cat(i, "\n")
  
  # Partitioning (clustering) of the data into k clusters “around medoids”
  clusteringList[[i-1]]   <- pam(distHarb, i)$clustering
  
  # Compute a number of distance based statistics, which can be used for cluster validation, 
  # comparison between clusterings and decision about the number of clusters
  clusterStatsList[[i-1]] <- fpc::cluster.stats(distHarb,clusteringList[[i-1]])
  #, G2 = T, G3 = T) #Too slow for these
}

avg.silwidths <- c(unlist(sapply(clusterStatsList, '[', 'avg.silwidth')))
widest.gaps   <- c(unlist(sapply(clusterStatsList, '[', 'widestgap')))
dunns         <- c(unlist(sapply(clusterStatsList, '[', 'dunn')))
dunn2s        <- c(unlist(sapply(clusterStatsList, '[', 'dunn2')))

png(
  "figures/harbison_choiceOfK_PAMcorrelations_silhouette.png",
  width = 500,
  height = 500
)
plot(
  2:20,
  avg.silwidths,
  type = "b",
  xlab = "Clusters",
  ylab = "",
  cex = 1.2,
  cex.lab = 1.5,
  cex.axis = 1.2
)
dev.off()

png(
  "figures/harbison_choiceOfK_PAMcorrelations_widestGap.png",
  width = 500,
  height = 500
)
plot(
  2:20,
  widest.gaps,
  type = "b",
  xlab = "Clusters",
  ylab = "",
  cex = 1.2,
  cex.lab = 1.5,
  cex.axis = 1.2
)
dev.off()

png(
  "figures/harbison_choiceOfK_PAMcorrelations_DunnsIndex.png",
  width = 500,
  height = 500
)
plot(
  2:20,
  dunns,
  type = "b",
  xlab = "Clusters",
  ylab = "",
  cex = 1.2,
  cex.lab = 1.5,
  cex.axis = 1.2
)
dev.off()

png(
  "figures/harbison_choiceOfK_PAMcorrelations_DunnsModifiedIndex.png",
  width = 500,
  height = 500
)
plot(
  2:20,
  dunn2s,
  type = "b",
  xlab = "Clusters",
  ylab = "",
  cex = 1.2,
  cex.lab = 1.5,
  cex.axis = 1.2
)
dev.off()

annotation_col = data.frame(K2 = factor(clusteringList[[1]]),
                            K3 = factor(clusteringList[[2]]),
                            K4 = factor(clusteringList[[3]]),
                            K5 = factor(clusteringList[[4]]),
                            K6 = factor(clusteringList[[5]]),
                            K7 = factor(clusteringList[[6]]),
                            K8 = factor(clusteringList[[7]]))

# On the basis of the diagnostics and the heatmap, it appears that 4 clusters
# is a reasonable number.

# Let's plot the data, showing the 8 cluster solution:
sortedClusters <- sort(clusteringList[[7]], index.return = T)

# Plot as a heatmap:
rowGaps <- c(match(2, sortedClusters$x), match(3, sortedClusters$x),
             match(4, sortedClusters$x), match(5, sortedClusters$x),
             match(6, sortedClusters$x), match(7, sortedClusters$x),
             match(8, sortedClusters$x)) - 1

png("figures/harbison_initialClusters_PAMgowers.png", height = 600, width = 800)
pheatmap(
    harbisonDataNumeric[sortedClusters$ix, ],
    cluster_rows = F,
    cluster_cols = T,
    gaps_row = rowGaps,
    fontsize_row = 5,
    cellheight = 0.5,
    show_rownames = F
)
dev.off()

# It seems like 8 clusters is a reasonable number
chosenK <- 8

pItem <- 0.95
n_reps <- 200
N <- dim(harbisonDataNumeric)[1]

coClusteringMatrix <- matrix(0, N, N)
rownames(coClusteringMatrix) <- colnames(coClusteringMatrix) <- rownames(harbisonDataNumeric)

normalisingMatrix <- matrix(0, N, N)
rownames(normalisingMatrix) <- colnames(normalisingMatrix) <- rownames(harbisonDataNumeric)

for(i in 1:n_reps){
  cat(i, "\n")
  set.seed(i)
  
  # Sample 80% of items
  items <- sample(N, floor(N*pItem))
  
  # Update counts of selected pairs of items
  normalisingMatrix <- normalisingMatrix + crossprod(t(c(1:N)%in%items))
  
  # Generate distances for selected items
  distHarbSubset <- as.dist(as.matrix(distHarb)[items, items])
  
  # Use PAM to cluster sampled items
  clusters <- pam(distHarbSubset, k = 8)$clustering
  
  for(j in 1:5){
    # Update counts of number of times each pair has been put in the same cluster
    coClusteringMatrix[items,items] <- coClusteringMatrix[items,items] +
      crossprod(t(as.numeric(clusters==j)))
  }
  
}

consensusMatrix <- coClusteringMatrix / normalisingMatrix

HharbisonData <- Heatmap(as.matrix(harbisonDataNumeric))
HconsensusMatrix <- Heatmap(consensusMatrix)
HconsensusMatrix + HharbisonData 

harbClusters <- clusteringList[[chosenK-1]]

save(harbClusters, consensusMatrix, file = "data/harbisonClusters_PAMgowers.RData")
write.table(
  as.data.frame(harbClusters),
  paste("goto-scores/harbison_", chosenK, "clusters_PAMgowers.csv", sep = ""),
  col.names = FALSE,
  quote = TRUE
)

#################################### Clustering all observations ###################################

# 4 clusters
harbison4Clusters <- clusteringList[[3]]
write.table(
  harbison4Clusters,
  "goto-scores/harbison_4clusters_PAMgowers.csv",
  col.names = FALSE,
  quote = TRUE
)

# 10 clusters
harbison10Clusters <- clusteringList[[9]]
write.table(
  harbison10Clusters,
  "goto-scores/harbison_10clusters_PAMgowers.csv",
  col.names = FALSE,
  quote = TRUE
)

# save(harbison4Clusters, harbison10Clusters, file = "harbisonClusters-PAMgowers.RData")

######################################### Plot clusters ############################################

inds <- match(rownames(harbisonDataNumeric), names(harbClusters))
harbClusters <- harbClusters[inds]

# Let's just plot the clusters:
table(harbClusters)
sortedClusters <- sort(harbClusters, index.return = T)

harbisonDataNumeric <-
  as.matrix(apply(harbisonDataNumeric, 2, as.factor))

col_fun = c("white", "black")

Hdata <-
  Heatmap(
    harbisonDataNumeric[sortedClusters$ix, ],
    cluster_rows = F,
    cluster_columns = T,
    row_split = sortedClusters$x,
    show_row_names = FALSE,
    show_column_names = FALSE,
    name = "Data",
    col = col_fun,
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
    )
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


# png("figures/harbison_initialClusters_PAMgowers_data.png",
#     height = 700,
#     width = 1000)
# Hdata + Hclusters
# dev.off()

png("figures/harbison_initialClusters_PAMgowers_consMatrix.png",
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
save(ccc, file = "cophenetic-correlation-coefficients/harbisonData_PAMgowers.RData")
  