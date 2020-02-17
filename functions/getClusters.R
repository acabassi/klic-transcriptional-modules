WhereToCut <- function(n) {
    attr(n, "height") <- 1
    if (!is.leaf(n)) {
        attr(n, "height") <- 2
        if (attr(n, "logEvidence") < 0) 
            attr(n, "height") <- 3
    }
    n
}

getClusters <- function (dendro, outputFile = "", verbose = FALSE) {
    
    # Apply WhereToCut to each node of a dendrogram recursively
    dendro <- dendrapply(dendro, WhereToCut)
    
    cutDendro <- cut(dendro, 2)
    
    
    nClusters <- length(cutDendro$lower)
    
    nTotalLabels <- length(labels(dendro))
    outputStrings  <- rep("", nTotalLabels)
    outputClusters <- vector(mode="numeric",length=nTotalLabels)
    counter <- 1
    for (i in 1:nClusters) {
        currentCluster <- cutDendro$lower[[i]]
        currentLabels <- labels(currentCluster)
        nLabels <- length(currentLabels)
        for (j in 1:nLabels) {
            outputStrings[counter] <- currentLabels[j]
            outputClusters[counter] <- i
            counter <- counter + 1
        }
    }
    names(outputClusters) <- outputStrings
    return(outputClusters)
}