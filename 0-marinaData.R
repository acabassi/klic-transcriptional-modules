################################# Transcriptional module discovery #################################
######################################### Expression data ##########################################

# We wish to cluster genes, to identify those that are functionally related.
# 
# We consider the example given in Section 3.3 of *Kirk, P., Griffin, J. E., Savage, R. S., 
# Ghahramani, Z., & Wild, D. L. (2012). Bayesian correlated clustering to integrate multiple 
# datasets. Bioinformatics (Oxford, England), 28(24), 3290–3297.*

rm(list=ls())

set.seed(1)

library(cluster)
library(pheatmap)
library(pracma)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
BiocManager::install(c("BHC", "ComplexHeatmap", "ConsensusClusterPlus"))
library(BHC)
library(ComplexHeatmap)
library(ConsensusClusterPlus)

source("functions/getClusters.R")

########################################### Load the data ##########################################

marinaDataNumeric <- read.table("data/marina_ppi_harbison.csv", 
  sep = ",", header = T, row.names = 1)

marinaData <- marinaDataNumeric

save(marinaData, file = "data/marinaData.RData")

#################################### Partitioning Around Medoids ###################################

cat("Marina data: Partitioning Around Medoids \n")

corMat <- cor(t(marinaData))
dissimilarity <- 1 - corMat

distMarina <- as.dist(dissimilarity)
clusteringList <- clusterStatsList <- list(length = 20)

for(i in 2:20){
  cat(i, "\n")
  
  # Partitioning (clustering) of the data into k clusters “around medoids”
  clusteringList[[i-1]]   <- pam(distMarina, i)$clustering
  
  # Compute a number of distance based statistics, which can be used for cluster validation, 
  # comparison between clusterings and decision about the number of clusters
  clusterStatsList[[i-1]] <- fpc::cluster.stats(distMarina,clusteringList[[i-1]])
                                                #, G2 = T, G3 = T) #Too slow for these
}

avg.silwidths <- c(unlist(sapply(clusterStatsList, '[', 'avg.silwidth')))
widest.gaps   <- c(unlist(sapply(clusterStatsList, '[', 'widestgap')))
dunns         <- c(unlist(sapply(clusterStatsList, '[', 'dunn')))
dunn2s        <- c(unlist(sapply(clusterStatsList, '[', 'dunn2')))

png(
  "figures/marina_choiceOfK_PAMcorrelations_silhouette.png",
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
  "figures/marina_choiceOfK_PAMcorrelations_widestGap.png",
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
  "figures/marina_choiceOfK_PAMcorrelations_DunnsIndex.png",
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
  "figures/marina_choiceOfK_PAMcorrelations_DunnsModifiedIndex.png",
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
                            K7 = factor(clusteringList[[6]]))

rownames(annotation_col) = rownames(corMat)

png("figures/marinaData_choiceOfK_heatmap.png", width = 700, height = 1000)
pheatmap(corMat, annotation_col = annotation_col, fontsize_row = 2,
fontsize_col = 2, cellheight = 1, cellwidth = 1)
dev.off()

# On the basis of the diagnostics and the heatmap, it appears that 4 clusters
# is a reasonable number.

# Let's plot the data, showing the 4 cluster solution:
sortedClusters <- sort(clusteringList[[3]], index.return = T)

# First plot as a heatmap:
rowGaps <- c(match(2, sortedClusters$x),
             match(3, sortedClusters$x),
             match(4, sortedClusters$x)) - 1
png("figures/marinaData_4Clusters_data.png", width = 700, height = 1000)
pheatmap(marinaDataNumeric[sortedClusters$ix,], cluster_rows = F, cluster_cols = T,
gaps_row = rowGaps, fontsize_row = 5, cellheight = 1)
dev.off()

# Next plot as time series:
png("figures/marina_4clusters_PAMcorrelations.png")
par(mfrow = c(2,2))
matplot(t(marinaDataNumeric[which(clusteringList[[3]]==1),]),
        type = "l", ylab = "", main = "Cluster 1")
matplot(t(marinaDataNumeric[which(clusteringList[[3]]==2),]),
        type = "l", ylab = "", main = "Cluster 2")
matplot(t(marinaDataNumeric[which(clusteringList[[3]]==3),]),
        type = "l", ylab = "", main = "Cluster 3")
matplot(t(marinaDataNumeric[which(clusteringList[[3]]==4),]),
        type = "l", ylab = "", main = "Cluster 4")
dev.off()

# It seems like 4 clusters is a reasonable number

pItem <- 0.8
n_reps <- 100
N <- dim(marinaData)[1]

consensusMatrix <- matrix(0, N, N)
rownames(consensusMatrix) <- colnames(consensusMatrix) <- rownames(marinaData)

normalisingMatrix <- matrix(0, N, N)
rownames(normalisingMatrix) <- colnames(normalisingMatrix) <- rownames(marinaData)

for(i in 1:n_reps){
  cat(i, "\n")
  set.seed(i)
  
  # Sample 80% of items
  items <- sample(N, floor(N*pItem))
  
  # Update counts of selected pairs of items
  normalisingMatrix <- normalisingMatrix + crossprod(t(c(1:N)%in%items))

  # Generate distances for selected items
  distMarina <- as.dist(dissimilarity[items, items])
  
  # Use PAM to cluster sampled items
  clusters <- pam(distMarina, k = 4)$clustering

  for(j in 1:4){
    # Update counts of number of times each pair has been put in the same cluster
    consensusMatrix[items,items] <- consensusMatrix[items,items] +
      crossprod(t(as.numeric(clusters==j)))
  }

}

consensusMatrix <- consensusMatrix / normalisingMatrix

HmarinaData <- Heatmap(as.matrix(marinaData))
HconsensusMatrix <- Heatmap(consensusMatrix)
HconsensusMatrix + HmarinaData 

marinaClusters <- clusteringList[[3]]

save(marinaClusters, consensusMatrix, file = "data/marinaClusters_PAMcorrelation.RData")
write.table(
    as.data.frame(marinaClusters),
    paste("goto-scores/marina_", 4, "clusters_PAMcorrelation.csv", sep = ""),
    col.names = FALSE,
    quote = TRUE
)

# Let's save the 10-cluster solution as well for the comparison with the other clustering algorithms

marina10Clusters <- clusteringList[[9]]
write.table(
    as.data.frame(marina10Clusters),
    paste("goto-scores/marina_", 10, "clusters_PAMcorrelation.csv", sep = ""),
    col.names = FALSE,
    quote = TRUE
)

######################################### Plot clusters ############################################

marClusters <- marinaClusters

inds <- match(rownames(marinaData), names(marClusters))
marClusters <- marClusters[inds]

# Let's just plot the clusters:
table(marClusters)
sortedClusters <- sort(marClusters, index.return = T)

col_fun = c("white", "black")


negative_numbers <- linspace(min(marinaDataNumeric), 0, n = 16)
positive_numbers <-
  linspace(0, max(marinaDataNumeric), n = ceil((
    max(marinaDataNumeric) / (negative_numbers[2] - negative_numbers[1])
  )))

my_colours <-
  colorRampPalette(c("#FF9900", "white"))(length(negative_numbers))
my_colours <-
  c(my_colours, colorRampPalette(c("white", "#146EB4"))(length(positive_numbers))[2:(length(positive_numbers))])

Hdata <-
  Heatmap(
    as.matrix(marinaData[sortedClusters$ix,]),
    col = my_colours,
    cluster_rows = F,
    cluster_columns = T,
    row_split = sortedClusters$x,
    show_row_names = FALSE,
    show_column_names = FALSE,
    name = "Data",
    width = unit(5, "cm"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 28),
      labels_gp = gpar(fontsize = 28),
      legend_direction = "horizontal"
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


# png(
#   "figures/marina_initialClusters_PAMcorrelations_data.png",
#   height = 700,
#   width = 1000
# )
# Hdata + Hclusters
# dev.off()

png(
  "figures/marina_initialClusters_PAMcorrelations_consMatrix.png",
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
save(ccc, file = "cophenetic-correlation-coefficients/marinaData_PAMcorrelation.RData")

############################################### Bonus ##############################################
################################# Bayesian Hierarchical Clustering #################################

# cat("Marina data: Bayesian Hierarchical Clustering \n")
# 
# bhc_marina20 <- list()
# time20 <- list()
# marinaClusters20 <- list()

### This takes a long time to run!
# 
# percentiles  <- FindOptimalBinning(marinaData, rownames(marinaData),
#                                    transposeData=TRUE, verbose=TRUE)
# discreteData <- DiscretiseData(t(marinaData), percentiles=percentiles)
# discreteData <- t(discreteData)
# 
# for(i in 1:25){
#   
#   cat(i, "\n")
#   set.seed(i)
#   
#   time20[[i]] <- system.time(bhc_marina20[[i]] <- bhc(discreteData, 
#                                             itemLabels = rownames(marinaData), 
#                                             randomised = T, m = 20, verbose = T,  
#                                             dataType = "time-course"))
#   marinaClusters20[[i]] <- getClusters(bhc_marina20[[i]])
# 
# save(list=ls(), file = "marinaClusters.RData")
# 
# }

# load("data/marinaAllClusters.RData") 

####################################### Co-clustering matrix #######################################

# n_reps <- 25
# N <- dim(marinaData)[1]
# 
# coClusteringMatrix <- matrix(0, N, N)
# rownames(coClusteringMatrix) <- colnames(coClusteringMatrix) <- rownames(marinaData)
# 
# for(i in 1:n_reps){
#   
#   if(i > 1){
#     differences <- sum(1-(marinaClusters20[[i]]-marinaClusters20[[i-1]]))
#     cat(differences, "\n")
#   } 
#   
#   maxK <- length(table(marinaClusters20[[i]]))
#   for(k in 1:maxK){
#     # Update counts of number of times each pair has been put in the same cluster
#     coClusteringMatrix <- coClusteringMatrix + crossprod(t(as.numeric(marinaClusters20[[i]]==k)))
#   }  
# }
# 
# coClusteringMatrix <- coClusteringMatrix/n_reps

########################################## Plot clusters ###########################################

# inds <- match(rownames(marinaData), names(marinaClusters20[[i]]))
# marinaClusters <- marinaClusters20[[i]][inds]
# 
# save(marinaClusters, coClusteringMatrix, file = "data/marinaClusters.RData")
# 
# # Let's just plot the clusters:
# table(marinaClusters)
# sortedClusters <- sort(marinaClusters, index.return = T)
# myBlues <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100)
# 
# # Plot as a heatmap:
# rowGaps <- c(match(2, sortedClusters$x), match(3, sortedClusters$x),
#              match(4, sortedClusters$x), match(5, sortedClusters$x)) - 1
# 
# col_fun = c("white", "black")
# 
# marinaData <- as.matrix(marinaData)
# 
# Hdata <- Heatmap(marinaData[sortedClusters$ix,], cluster_rows = F, 
#                  cluster_columns = T, row_split = sortedClusters$x, 
#                  show_row_names = FALSE, name = "Data")
# HconsMatrix <- Heatmap(coClusteringMatrix[sortedClusters$ix,sortedClusters$ix], 
#                        # cluster_columns = FALSE,  
#                        show_column_names = FALSE, name = "Similarities", col = myBlues, 
#                        row_split = sortedClusters$x, 
#                        column_split = sortedClusters$x)
# Hclusters <- Heatmap(as.matrix(as.factor(sortedClusters$x)), show_row_names = FALSE, 
#                      name = "Clusters")
# 
# # Hdata + HconsMatrix + Hclusters # Too big
# 
# pdf("figures/marinaDataClusters.pdf")
# Hdata + Hclusters
# dev.off()
# 
# pdf("figures/marinaConsMatrixClusters")
# HconsMatrix + Hclusters
# dev.off()
