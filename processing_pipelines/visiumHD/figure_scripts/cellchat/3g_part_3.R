#remotes::install_github("jinworks/CellChat@534b91d")

suppressPackageStartupMessages({
    library(CellChat)
    library(zellkonverter)
    library(SingleCellExperiment())
})

cellchat <- readRDS("tmp/cellchat.rds")

# change the order of the labels
cellchat@meta$labels <- forcats::fct_relevel(cellchat@meta$labels, "Cd8 top", "Cd8 crypt", "Cd8 muscularis")
cellchat <- updateClusterLabels(cellchat, new.order = levels(cellchat@meta$labels))

#pathways.show <- c("TGFb")
#netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

library(NMF)
#selectK(cellchat, pattern = "outgoing")

nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

library(ggalluvial)
netAnalysis_river(cellchat, pattern = "incoming")

#netAnalysis_dot(cellchat, pattern = "incoming", group.show = c("Cd8 top", "Cd8 crypt", "Cd8 muscularis"), font.size = 10)

#write to csv
write.csv(cellchat@netP$pattern$incoming$data, file = "/tmp/incoming.csv")
#write to csv
write.csv(cellchat@netP$pattern$outgoing$data, file = "/tmp/outgoing.csv")

signaling_pathway <- c('CXCL')
netAnalysis_contribution(cellchat, signaling = signaling_pathway)



# Get all available signaling pathways
signaling_pathways <- cellchat@netP$pathways

# Initialize an empty data frame to store contributions from all pathways
df_all <- data.frame()

# Loop over all signaling pathways
for (signaling_pathway in signaling_pathways) {

  # Perform the analysis for the current pathway
  pairLR <- searchPair(signaling = signaling_pathway, pairLR.use = cellchat@LR$LRsig, key = "pathway_name", matching.exact = TRUE, pair.only = TRUE)

  # Extract interaction names for the current pathway
  pair.name.use <- select(cellchat@DB$interaction[rownames(pairLR), ], "interaction_name_2")

  # Ensure pair names are extracted correctly
  net <- cellchat@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  print(pairLR.name)
  # Subset prob and pval matrices to include only interactions from the current pathway
  prob <- net$prob[, , pairLR.name]
  pval <- net$pval[, , pairLR.name]

  # Threshold for significance
  thresh <- 0.05
  prob[pval > thresh] <- 0

  # Handle the case where there's only one member in the pathway
  if (length(dim(prob)) == 3) {
    # If prob is 3D, apply sum across the third dimension
    pSum <- apply(prob, 3, sum)
  } else if (length(dim(prob)) == 2) {
    # If prob is 2D, just sum the matrix
    pSum <- sum(prob)
    pSum <- c(pSum)  # Convert to a vector to match the structure
  } else {
    # In case of unexpected dimensions, set pSum to 0
    pSum <- 0
  }

  pSum.max <- sum(prob)
  pSum <- pSum / pSum.max
  pSum[is.na(pSum)] <- 0

  # Ensure pair names are mapped correctly
  if (!is.null(pairLR.use.name)) {
    pair.name <- pair.name.use[, "interaction_name_2"]
    pair.name <- factor(pair.name, levels = unique(pair.name))
  }

  #print(pair.name.use)
  # Create a data frame for the current pathway's contribution
  df_current <- data.frame(name = signaling_pathway, contribution = pSum)
  df_current$interaction <- pair.name

  # Handle missing values
  df_current$contribution[is.na(df_current$contribution)] <- 0

  # Combine current pathway results with the overall data frame
  df_all <- rbind(df_all, df_current)
}

# Ensure the names are ordered by contribution
df_all <- df_all[order(df_all$contribution, decreasing = TRUE), ]


# Save the df_all dataframe as a CSV file
write.csv(df_all, file = "/figures/visiumHD_pathway_contributions.csv", row.names = F)
