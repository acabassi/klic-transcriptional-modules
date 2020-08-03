################################ Transcriptional module discovery ##################################
############################################# KLIC #################################################

rm(list = ls())

set.seed(1)

library(ComplexHeatmap)
library(colorspace)
library(rcartocolor)

# library(devtools)
# install_github("acabassi/coca")
library(coca)
# install_github("acabassi/klic")
library(klic)

########################################## Load ChIP data ##########################################

harbisonClusteringAlgorithm <- "PAMgowers" 
# Must be one of "GBNP", "BHC", "PAMgowers".

load(paste(
    "data/harbisonClusters_",
    harbisonClusteringAlgorithm,
    ".RData",
    sep = ""
))
ccH <- consensusMatrix
rm(consensusMatrix)

######################################### Load Marina data #########################################
load("data/marinaClusters_PAMcorrelation.RData")

sum(1 - (rownames(ccH) == rownames(consensusMatrix)))
sum(1 - (rownames(consensusMatrix) == rownames(ccH)))

# -------- Only run this if the previous line gives "5" as output ------------
# odd_names <- names(marinaClusters)[grep("-A", names(marinaClusters))]
# rownames(ccH)[which(!rownames(ccH)==rownames(consensusMatrix))] <- odd_names
# colnames(ccH)[which(!colnames(ccH)==colnames(consensusMatrix))] <- odd_names
# ----------------------------------------------------------------------------

ccM <- consensusMatrix
rm(consensusMatrix)

########################################## Spectrum shift ##########################################

# Cophenetic correlation before spectrum shift
cophCorr_ccH_before <- copheneticCorrelation(ccH)
cophCorr_ccM_before <- copheneticCorrelation(ccM)

ccH <- spectrumShift(ccH,    verbose = TRUE, coeff = 1.01)
ccM <- spectrumShift(ccM, verbose = TRUE, coeff = 1.01)

cophCorr_ccH_after <- copheneticCorrelation(ccH)
cophCorr_ccM_after <- copheneticCorrelation(ccM)

###################################### Create array of all matrices ################################

n_datasets <- 2
datasetNames <- c("Harbison", "Marina")
N <- dim(ccH)[1]
CM <- array(NA, c(N, N, n_datasets))
dimnames(CM) <- list(rownames(ccH), rownames(ccH), datasetNames)

CM[rownames(ccH), colnames(ccH), "Harbison"] <- ccH
CM[rownames(ccM), colnames(ccM), "Marina"]   <- ccM

################################################# KLIC #############################################

# Use localised multiple kernel k-means to integrate the datasets

# Initialise array of kernel matrices

maxK <- 25
KM <- array(0, c(N, N, maxK - 1))
clLabels <- array(NA, c(maxK - 1, N))
theta <- array(NA, c(maxK - 1, N, n_datasets))

parameters <- list()
parameters$iteration_count <-
    100 # set the maximum number of iterations

for(i in 2:maxK) {
    # Use kernel k-means with K=i to find weights and cluster labels
    parameters$cluster_count <- i # set the number of clusters K
    lmkkm <- lmkkmeans_missingData(CM, parameters, verbose = TRUE)

    # Compute weighted matrix
    for (j in 1:dim(CM)[3]) {
        KM[, , i - 1] <-
            KM[, , i - 1] + (lmkkm$Theta[, j] %*% t(lmkkm$Theta[, j])) * CM[, , j]
    }

    # Save cluster labels
    clLabels[i - 1, ] <- lmkkm$clustering
    theta[i - 1, , ] <- lmkkm$Theta
}

######################################### Maximise silhouette ######################################

maxSil <- maximiseSilhouette(KM, clLabels, maxK = maxK)

save(
    KM,
    clLabels,
    theta,
    maxSil,
    file = paste(
        "data/klic_output_",
        harbisonClusteringAlgorithm,
        ".RData",
        sep = ""
    )
)

load(paste(
    "data/klic_output_",
    harbisonClusteringAlgorithm,
    ".RData",
    sep = ""
))

dimnames(theta) <- list(paste(as.character(2:maxK), "_clusters"),
                        rownames(CM),
                        c("ChIP", "Expression")
)


############################################ Plot output ###########################################

bestK <- maxSil$K
inds <- rownames(CM)
clustersBestK <- clLabels[bestK - 1, ]
names(clustersBestK) <- rownames(CM)
# names(clustersBestK)[which(names(clustersBestK)%in%replacement_names)] <- odd_names

# Write clusters to csv files:
write.table(
    as.data.frame(clustersBestK),
    paste(
        "goto-scores/klic_",
        bestK,
        "clusters_",
        harbisonClusteringAlgorithm,
        ".csv",
        sep = ""
    ),
    col.names = FALSE,
    quote = TRUE
)

# This is just to make sure that our GOTO scores are better than what we would
# have obtained just by chance 
# shuffledClustersBestK <- sample(clustersBestK)
# table(clustersBestK)
# table(shuffledClustersBestK)
# names(shuffledClustersBestK) <- names(clustersBestK)
# sum(1-(clustersBestK==shuffledClustersBestK))
# write.table(
#     as.data.frame(shuffledClustersBestK),
#     paste(
#         "goto-scores/klic_",
#         bestK,
#         "clusters_shuffled_",
#         harbisonClusteringAlgorithm,
#         ".csv",
#         sep = ""
#     ),
#     col.names = FALSE,
#     quote = TRUE
# )

# Let's just plot the clusters:
table(clustersBestK)
sortedClusters <- sort(clustersBestK, index.return = T)
myBlues <-
    colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)

# Color blind friendly palette for the annotations
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

my_anno_legend_param <- function(name, nrow=1) {
    return(list(
        title = name,
        labels_gp = gpar(fontsize = 18),
        title_gp = gpar(fontsize = 22),
        nrow = nrow
    ))
}

my_heatmap_legend_param <- list(labels_gp = gpar(fontsize = 18),
                                title_gp = gpar(fontsize = 22),
                                direction = "horizontal",
                                nrow = 1,
                                legend_width = unit(4.5, "cm"))
col_bw = c("white", "black")

load("data/harbisonData.RData")
harbCl <- as.factor(as.character(harbClusters)) 
names(harbCl) <- names(harbClusters)
# inds <- match(rownames(harbisonDataNumeric), names(harbCl))
# harbCl <- harbCl[inds]

nColor <- 12
bigPalette <- carto_pal(nColor, "Safe")

HharClusters <- rowAnnotation(
    harbisonClusters = harbCl,
    name = "Harbison Clusters",
    annotation_legend_param = my_anno_legend_param("ChIP clusters", nrow = 2),
    show_annotation_name = FALSE,
    col = list(
        harbisonClusters = c(
            "1" = bigPalette[1],
            "2" = bigPalette[2],
            "3" = bigPalette[3],
            "4" = bigPalette[4],
            "5" = bigPalette[5],
            "6" = bigPalette[6],
            "7" = bigPalette[7],
            "8" = bigPalette[8],
            "9" = bigPalette[9],
            "10" = bigPalette[10],
            "11" = bigPalette[11],
            "13" = bigPalette[12] # "GBNP", cluster 12 is empty
        )
    )
)

HharbisonData <-
    Heatmap(
        as.matrix(harbisonDataNumeric),
        cluster_rows = F,
        cluster_columns = T,
        row_split = sortedClusters$x,
        row_gap = unit(4, "mm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "ChIP data",
        right_annotation = HharClusters,
        col = col_bw,
        heatmap_legend_param = my_heatmap_legend_param,
        width = unit(8, "cm")
    )
HharbisonData

load("data/marinaData.RData")
HmarClusters <-
    rowAnnotation(
        MarinaClusters = as.factor(marinaClusters),
        name = "Exp. clusters",
        annotation_legend_param = my_anno_legend_param("Exp. clusters"),
        show_annotation_name = FALSE,
        col = list(
            MarinaClusters = c(
                "1" = bigPalette[1],
                "2" = bigPalette[2],
                "3" = bigPalette[3],
                "4" = bigPalette[4],
                "5" = bigPalette[5],
                "6" = bigPalette[6],
                "7" = bigPalette[7],
                "8" = bigPalette[8],
                "9" = bigPalette[9],
                "10" = bigPalette[10],
                "11" = bigPalette[11],
                "12" = bigPalette[12]
            )
        )
    )
HmarinaData <- Heatmap(
    as.matrix(marinaData),
    cluster_rows = F,
    cluster_columns = T,
    row_split = sortedClusters$x,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    name = "Exp. data",
    right_annotation = HmarClusters,
    heatmap_legend_param = my_heatmap_legend_param,
    width = unit(4.1, "cm")
)
HmarinaData

Hclusters <-
    rowAnnotation(
        finalClusters = as.factor(sortedClusters$x),
        name = "Final clusters",
        annotation_legend_param = my_anno_legend_param("KLIC clusters"),
        show_annotation_name = FALSE,
        col = list(
            finalClusters = c(
                "1" = cbPalette[1],
                "2" = cbPalette[2],
                "3" = cbPalette[3],
                "4" = cbPalette[4],
                "5" = cbPalette[5]
            )
        )
    )
HwMatrix <-
    Heatmap(
        KM[sortedClusters$ix, sortedClusters$ix, bestK - 1],
        # cluster_columns = FALSE,
        show_column_names = FALSE,
        name = "Weighted kernel",
        col = myBlues,
        row_split = sortedClusters$x,
        column_split = sortedClusters$x,
        right_annotation = Hclusters,
        heatmap_legend_param = my_heatmap_legend_param,
        width = unit(20, "cm"),
        height = unit(20, "cm")
    )

png(
    paste0(
        "figures/klic_weightedKernel_", bestK, "clusters_",
        harbisonClusteringAlgorithm,
        ".png"
    ),
    height = 800,
    width = 900
)
draw(
    HwMatrix,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

png(
    paste(
        "figures/klic_harbisonData_finalClusters_",
        harbisonClusteringAlgorithm,
        ".png",
        sep = ""
    ),
    height = 800, 
    width = 1200
)
HwMatrix + HharbisonData
dev.off()

png(
    paste(
        "figures/klic_marinaData_finalClusters_",
        harbisonClusteringAlgorithm,
        ".png",
        sep = ""
    ),
    height = 800,
    width = 1200,
)
HwMatrix + HmarinaData
dev.off()

png(
    paste(
        "figures/klic_finalClusters_",
        harbisonClusteringAlgorithm,
        ".png",
        sep = ""
    ),
    height = 800,
    width = 1400
)
draw(
    HwMatrix + HharbisonData + HmarinaData,
    merge_legend = TRUE,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

HmarinaData + HharbisonData

# Plot weights

weights_palette <- rev(sequential_hcl(10, palette = "Mint"))
Hweights <-
    Heatmap(
        theta[bestK, , ],
        col = weights_palette,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        # right_annotation = Hclusters,
        name = "Weights",
        heatmap_legend_param = my_heatmap_legend_param,
        column_names_side = "top",
        # column_names_rot = 0,
        column_names_gp = gpar(fontsize = 18),
        # column_names_centered = TRUE,
        width = unit(1, "cm"),
    )

png(
    paste(
        "figures/klic_weights_",
        harbisonClusteringAlgorithm,
        ".png",
        sep = ""
    ),
    height = 800,
    width = 1400
)
draw(
    HwMatrix + Hweights + HharbisonData + HmarinaData ,
    merge_legend = TRUE,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

############## Repeat the analysis with K = 10, for comparison with the other methods ##############

bestK <- 10

clustersBestK <- clLabels[bestK - 1, ]
names(clustersBestK) <- rownames(CM)
# names(clustersBestK)[which(names(clustersBestK)%in%replacement_names)] <- odd_names

# Write clusters to csv files:
write.table(
    as.data.frame(clustersBestK),
    paste(
        "goto-scores/klic_",
        bestK,
        "clusters_",
        harbisonClusteringAlgorithm,
        ".csv",
        sep = ""
    ),
    col.names = FALSE,
    quote = TRUE
)

shuffledClustersBestK <- sample(clustersBestK)
names(shuffledClustersBestK) <- names(clustersBestK)
table(clustersBestK)
table(shuffledClustersBestK)
sum(1-(clustersBestK==shuffledClustersBestK))
write.table(
    as.data.frame(shuffledClustersBestK),
    paste("goto-scores/klic_", bestK, "clusters_shuffled.csv", sep = ""),
    col.names = FALSE,
    quote = TRUE
)

# Let's just plot the clusters:
sortedClusters <- sort(clustersBestK, index.return = T)

inds <- match(rownames(harbisonDataNumeric), names(harbCl))
harbCl <- harbCl[inds]

HharClusters <- rowAnnotation(
    harbisonClusters = harbCl,
    name = "Harbison Clusters",
    annotation_legend_param = my_anno_legend_param("Harbison Clusters"),
    show_annotation_name = FALSE
)

HharbisonData <-
    Heatmap(
        as.matrix(harbisonDataNumeric),
        cluster_rows = F,
        cluster_columns = T,
        row_split = sortedClusters$x,
        row_gap = unit(4, "mm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "Harbison Data",
        right_annotation = HharClusters,
        col = col_bw,
        heatmap_legend_param = my_heatmap_legend_param,
        width = unit(8, "cm")
    )
HharbisonData

load("data/marinaData.RData")
HmarClusters <-
    rowAnnotation(
        MarinaClusters = as.factor(marinaClusters),
        name = "Marina clusters",
        annotation_legend_param = my_anno_legend_param("Marina Clusters"),
        show_annotation_name = FALSE,
        col = list(
            MarinaClusters = c(
                "1" = cbPalette[5],
                "2" = cbPalette[6],
                "3" = cbPalette[7],
                "4" = cbPalette[8]
            )
        )
    )
HmarinaData <- Heatmap(
    as.matrix(marinaData),
    cluster_rows = F,
    cluster_columns = T,
    row_split = sortedClusters$x,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    name = "Marina Data",
    right_annotation = HmarClusters,
    heatmap_legend_param = my_heatmap_legend_param,
    width = unit(4.1, "cm")
)
HmarinaData

Hclusters <-
    rowAnnotation(
        finalClusters = as.factor(sortedClusters$x),
        name = "Final clusters",
        annotation_legend_param = my_anno_legend_param("KLIC Clusters"),
        show_annotation_name = FALSE
        )

HwMatrix <-
    Heatmap(
        KM[sortedClusters$ix, sortedClusters$ix, bestK - 1],
        # cluster_columns = FALSE,
        show_column_names = FALSE,
        name = "Weighted kernel",
        col = myBlues,
        row_split = sortedClusters$x,
        column_split = sortedClusters$x,
        right_annotation = Hclusters,
        heatmap_legend_param = my_heatmap_legend_param
    )

# HwMatrix
# 
# HwMatrix + HharbisonData
# 
# HwMatrix + HmarinaData
# 
# HwMatrix + HharbisonData + HmarinaData
# 
# HmarinaData + HharbisonData

png(
    paste(
        "figures/klic_10Clusters_",
        harbisonClusteringAlgorithm,
        ".png",
        sep = ""
    ),
    height = 800,
    width = 1200
)
draw(
    HwMatrix + HharbisonData + HmarinaData,
    merge_legend = TRUE,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()
