################################ Transcriptional module discovery ##################################
############################################# COCA #################################################

rm(list = ls())

set.seed(1)

library(coca)

load("data/marinaClusters_PAMcorrelation.RData")
rm(consensusMatrix)

harbisonClusteringAlgorithm <- "PAMgowers" 
# Must be one of "GBNP", "BHC", "PAMgowers".

load(paste("data/harbisonClusters_", harbisonClusteringAlgorithm,".RData", sep=""))
rm(consensusMatrix)

n_datasets <- 2

sum(1 - names(harbClusters) %in% names(marinaClusters))
sum(1 - names(marinaClusters) %in% names(harbClusters))

# -------- Only run this if the previous line gives "5" as output ----
# odd_names <- names(marinaClusters)[grep("-A", names(marinaClusters))]
# remove_part = function(string){
#     new_string <- substr(string,1,nchar(string)-2)
#     return(new_string)
# }
# replacement_names <- unlist(lapply(odd_names, remove_part))
# names(marinaClusters)[grep("-A", names(marinaClusters))] <- replacement_names
# sum(1 - names(harbClusters) %in% names(marinaClusters))
# sum(1 - names(marinaClusters) %in% names(harbClusters))
# --------------------------------------------------------------------

moc <- matrix(0, length(marinaClusters), n_datasets)

rownames(moc) <- names(harbClusters)
colnames(moc) <- c("Harbison", "Marina") #, "PPI"
moc[names(harbClusters),   "Harbison"] <- harbClusters
moc[names(marinaClusters), "Marina"]   <- marinaClusters

expandedMOC <- expandMOC(moc, datasetNames = colnames(moc))
plotMOC(
    expandedMOC$moc,
    datasetIndicator = expandedMOC$datasetIndicator,
    datasetNames = colnames(moc)[expandedMOC$datasetIndicator]
)

moc <- expandedMOC$moc
coca <-
    coca(moc,
         maxK = 10,
         verbose = TRUE,
         returnAllMatrices = TRUE)

annotations <- as.data.frame(as.factor(coca$clusterLabels))
names(annotations) <- "Clusters"
rownames(annotations) <- rownames(moc)

dev.off()
plotMOC(
    moc,
    datasetIndicator = expandedMOC$datasetIndicator,
    datasetNames = colnames(moc)[expandedMOC$datasetIndicator],
    annotations = annotations,
    save = TRUE,
    fileName = paste(
        "figures/coca_MOC_10clusters_",
        harbisonClusteringAlgorithm,
        ".png",
        sep = ""
    )
)

# Save clusters to csv
# rownames(annotations)[which(rownames(annotations) %in% replacement_names)] <-
#     odd_names
annotations <- as.data.frame(annotations)
annotations$Clusters <- as.numeric(annotations$Clusters)
write.table(
    as.data.frame(annotations),
    paste(
        "goto-scores/coca_",
        coca$K,
        "clusters_",
        harbisonClusteringAlgorithm,
        ".csv",
        sep = ""
    ),
    col.names = FALSE,
    quote = TRUE
)

################# Repeat the same steps with K equal to value selected by KLIC  ####################

n_clusters <- 3
coca <- coca(moc,
             K = n_clusters,
             verbose = TRUE,
             returnAllMatrices = TRUE)

annotations <- as.data.frame(as.factor(coca$clusterLabels))
names(annotations) <- "Clusters"
rownames(annotations) <- rownames(moc)

dev.off()
plotMOC(
    moc,
    datasetIndicator = expandedMOC$datasetIndicator,
    datasetNames = colnames(moc)[expandedMOC$datasetIndicator],
    annotations = annotations,
    save = TRUE,
    fileName = paste(
        "figures/coca_MOC_",
        n_clusters,
        "clusters_",
        harbisonClusteringAlgorithm,
        ".png",
        sep = ""
    )
)

# Save clusters to csv
# rownames(annotations)[which(rownames(annotations) %in% replacement_names)] <-
    # odd_names
annotations <- as.data.frame(annotations)
annotations$Clusters <- as.numeric(annotations$Clusters)
write.table(
    as.data.frame(annotations),
    paste(
        "goto-scores/coca_",
        n_clusters,
        "clusters_",
        harbisonClusteringAlgorithm,
        ".csv",
        sep = ""
    ),
    col.names = FALSE,
    quote = TRUE
)
