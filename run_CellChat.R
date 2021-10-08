# load libraries
library(CellChat)
library(patchwork)
library(Hmisc)
library(igraph)

# set the working directory
data.dir <- './'
setwd(data.dir)


# load data
counts.ctrl.00d <- read.table("counts_control_00d.txt", header = T)
counts.ctrl.00d <- normalizeData(counts.ctrl.00d)
#rownames(counts.ctrl.00d) <- counts.ctrl.00d[,1]
#counts.ctrl.00d <- counts.ctrl.00d[,-1]
meta.ctrl.00d <- read.table("meta_control_00d.txt", header = T)
counts.ctrl.07d <- read.table("counts_control_07d.txt", header = T)
counts.ctrl.07d <- normalizeData(counts.ctrl.07d)
#rownames(counts.ctrl.07d) <- counts.ctrl.00d[,1]
#counts.ctrl.07d <- counts.ctrl.07d[,-1]
meta.ctrl.07d <- read.table("meta_control_07d.txt", header = T)
counts.ctrl.28d <- read.table("counts_control_28d.txt", header = T)
counts.ctrl.28d <- normalizeData(counts.ctrl.28d)
meta.ctrl.28d <- read.table("meta_control_28d.txt", header = T)
#rownames(counts.ctrl.28d) <- counts.ctrl.28d[,1]
#counts.ctrl.28d <- counts.ctrl.28d[,-1]

counts.plx.00d <- read.table("counts_PLX5622_00d.txt", header = T)
counts.plx.00d <- normalizeData(counts.plx.00d)
#rownames(counts.plx.00d) <- counts.plx.00d[,1]
#counts.plx.00d <- counts.plx.00d[,-1]
meta.plx.00d <- read.table("meta_PLX5622_00d.txt", header = T)
counts.plx.07d <- read.table("counts_PLX5622_07d.txt", header = T)
counts.plx.07d <- normalizeData(counts.plx.07d)
#rownames(counts.plx.07d) <- counts.plx.07d[,1]
#counts.plx.07d <- counts.plx.07d[,-1]
meta.plx.07d <- read.table("meta_PLX5622_07d.txt", header = T)
counts.plx.28d <- read.table("counts_PLX5622_28d.txt", header = T)
counts.plx.28d <- normalizeData(counts.plx.28d)
#rownames(counts.plx.28d) <- counts.plx.28d[,1]
#counts.plx.28d <- counts.plx.28d[,-1]
meta.plx.28d <- read.table("meta_PLX5622_28d.txt", header = T)


# the function to use CellChat
use_CellChat <- function(counts, meta)
  # counts : the gene expression matrix
  # meta   : the meta data matched with the gene expression matrix
{
  counts <- as.matrix(counts) # transform data frame into matrix
  cellchat <- createCellChat(object = counts, meta = meta, group.by = "Cell_type") # create a CellChat object
  levels(cellchat@idents) # check the cell types
  CellChatDB <- CellChatDB.mouse # set the database according to the organism
  dplyr::glimpse(CellChatDB$interaction) # show the structure of the database
  unique(CellChatDB$interaction$annotation) # check the annotation of interactions
  # CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
  CellChatDB.use <- CellChatDB # use all databases
  cellchat@DB <- CellChatDB.use # set the used database in the object
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- computeCommunProb(cellchat, raw.use = T, nboot = 1, population.size = T)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat) # get the pairs of cell types or ligand receptors in the interactions
  cell_type_pair <- paste0(df.net$source, '_', df.net$target) # build the pairs of interaction cell types
  mean_exp <- sapply(1:nrow(df.net), function(x) {
    ct1 <- df.net$source[x] # source cell type
    ct2 <- df.net$target[x] # target cell type
    cells1 <- meta[meta$Cell_type == ct1,][[1]] # cells belonging to the source cell type
    cells2 <- meta[meta$Cell_type == ct2,][[1]] # cells belonging to the target cell type
    # gene_pairs <- strsplit(df.net$interaction_name_2[x], split = ' - ', fixed = T)
    # genes1 <- strsplit(gene_pairs[[1]][1], split = '+', fixed = T)[[1]]
    # genes2 <- strsplit(gene_pairs[[1]][2], split = '+', fixed = T)[[1]]
    # genes1 <- gsub(pattern = '(', replacement = '', genes1, fixed = T)
    # genes1 <- gsub(pattern = ')', replacement = '', genes1, fixed = T)
    # genes1 <- gsub(pattern = ' ', replacement = '', genes1, fixed = T)
    # genes2 <- gsub(pattern = '(', replacement = '', genes2, fixed = T)
    # genes2 <- gsub(pattern = ')', replacement = '', genes2, fixed = T)
    # genes2 <- gsub(pattern = ' ', replacement = '', genes2, fixed = T)
    gene_pair1 <- df.net$ligand[x] # genes corresponding to the ligand
    genes1 <- capitalize(tolower(strsplit(gene_pair1, split = '_')[[1]]))
    # genes1 <- sapply(genes1, capitalize)
    gene_pair2 <- df.net$receptor[x] # genes corresponding to the receptor
    genes2 <- capitalize(tolower(strsplit(gene_pair2, split = '_')[[1]]))
    # genes2 <- sapply(genes2, capitalize)
    if (length(genes2) > 1 & genes2[2] == 'R2')
    {
      genes2[2] <- gsub(pattern = '1', replacement = '2', genes2[1])
    }
    cat("Ligand genes: ", genes1, "\n")
    cat("Receptor genes: ", genes2, "\n\n")
    mean1 <- mean(counts[genes1, cells1]) # the mean expression value of the ligand within the source cell type
    mean2 <- mean(counts[genes2, cells2]) # the mean expression value of the receptor within the source cell type
    avg <- (mean1 + mean2) / 2 # the mean expression value of the ligands and receptors
    return(avg)
  })
  df.net <- cbind(cell_type_pair, mean_exp, df.net) # add the mean expression and cell type pairs
  return(df.net)
  # cellchat <- computeCommunProbPathway(cellchat)
  # 
  # return(cellchat)
}


generateObject <- function(data, meta)
{
  data <- as.matrix(data)
  data <- as(data, "dgCMatrix")
  new.meta <- meta[,-1]
  rownames(new.meta) <- meta[, 1]
  colnames(new.meta) <- c("labels", "time", "condition")
  obj <- list(data = data, meta = new.meta)
  return(obj)
}


generate_dot_plot <- function(df.net, title)
  # df.net : the data frame to represent the interaction network
{
  p <- ggplot(df.net, aes(x = cell_type_pair, y = interaction_name_2, color = prob, size = mean_exp)) + 
    geom_point() + theme(axis.text.x = element_text(angle = 90)) + 
    xlab("Cell type pairs") + ylab("L-R pairs") + ggtitle(title)
  ggsave(p, filename = paste0(title, '.pdf'))
}


# visualization

# visualize the dot plots for the control and PLX data
# pdf("control_dot_plots.pdf")
# par(mfrow = c(2, 3), xpd=TRUE)
#layout(matrix(c(1, 1, 1, 2, 2, 1, 2, 2, 3, 1, 3, 2), 3, 2, byrow = TRUE))

# control 00d
net.ctrl.00d <- use_CellChat(counts.ctrl.00d, meta.ctrl.00d)
generate_dot_plot(net.ctrl.00d, "Control 00d")

# control 07d
net.ctrl.07d <- use_CellChat(counts.ctrl.07d, meta.ctrl.07d)
generate_dot_plot(net.ctrl.07d, "Control 07d")

# control 28d
net.ctrl.28d <- use_CellChat(counts.ctrl.28d, meta.ctrl.28d)
generate_dot_plot(net.ctrl.28d, "Control 28d")

# PLX 00d
net.plx.00d <- use_CellChat(counts.plx.00d, meta.plx.00d)
generate_dot_plot(net.plx.00d, "PLX 00d")

# PLX 07d
net.plx.07d <- use_CellChat(counts.plx.07d, meta.plx.07d)
generate_dot_plot(net.plx.07d, "PLX 07d")

# PLX 28d
net.plx.28d <- use_CellChat(counts.plx.28d, meta.plx.28d)
generate_dot_plot(net.plx.28d, "PLX 28d")

# dev.off()


# # visualize the aggregated cell-cell communication network for the control data
# pdf("control_aggregate_net.pdf")
# par(mfrow = c(3, 2), xpd=TRUE)
# #layout(matrix(c(1, 1, 1, 2, 2, 1, 2, 2, 3, 1, 3, 2), 3, 2, byrow = TRUE))
# 
# # control 00d
# # cellchat.ctrl.00d <- use_CellChat(counts.ctrl.00d, meta.ctrl.00d)
# net.ctrl.00d <- use_CellChat(counts.ctrl.00d, meta.ctrl.00d)
# cellchat.ctrl.00d <- aggregateNet(cellchat.ctrl.00d)
# groupSize.ctrl.00d <- as.numeric(table(cellchat.ctrl.00d@idents))
# netVisual_circle(cellchat.ctrl.00d@net$count, vertex.weight = groupSize.ctrl.00d, weight.scale = T, 
#                  label.edge= F, title.name = "Number of interactions - 00d")
# netVisual_circle(cellchat.ctrl.00d@net$weight, vertex.weight = groupSize.ctrl.00d, weight.scale = T, 
#                  label.edge= F, title.name = "Interaction weights/strength - 00d")
# 
# # control 07d
# cellchat.ctrl.07d <- use_CellChat(counts.ctrl.07d, meta.ctrl.07d)
# cellchat.ctrl.07d <- aggregateNet(cellchat.ctrl.07d)
# groupSize.ctrl.07d <- as.numeric(table(cellchat.ctrl.07d@idents))
# netVisual_circle(cellchat.ctrl.07d@net$count, vertex.weight = groupSize.ctrl.07d, weight.scale = T, 
#                  label.edge= F, title.name = "Number of interactions - 07d")
# netVisual_circle(cellchat.ctrl.07d@net$weight, vertex.weight = groupSize.ctrl.07d, weight.scale = T, 
#                  label.edge= F, title.name = "Interaction weights/strength - 07d")
# 
# # control 28d
# cellchat.ctrl.28d <- use_CellChat(counts.ctrl.28d, meta.ctrl.28d)
# cellchat.ctrl.28d <- aggregateNet(cellchat.ctrl.28d)
# groupSize.ctrl.28d <- as.numeric(table(cellchat.ctrl.28d@idents))
# netVisual_circle(cellchat.ctrl.28d@net$count, vertex.weight = groupSize.ctrl.28d, weight.scale = T, 
#                  label.edge= F, title.name = "Number of interactions - 28d")
# netVisual_circle(cellchat.ctrl.28d@net$weight, vertex.weight = groupSize.ctrl.28d, weight.scale = T, 
#                  label.edge= F, title.name = "Interaction weights/strength - 28d")
# 
# dev.off()
# 
# # visualize the aggregated cell-cell communication network for the PLX5622 data
# pdf("PLX5622_aggregate_net.pdf")
# par(mfrow = c(3, 2), xpd=TRUE)
# #layout(matrix(c(1, 1, 1, 2, 2, 1, 2, 2, 3, 1, 3, 2), 3, 2, byrow = TRUE))
# 
# # control 00d
# cellchat.plx.00d <- use_CellChat(counts.plx.00d, meta.plx.00d)
# cellchat.plx.00d <- aggregateNet(cellchat.plx.00d)
# groupSize.plx.00d <- as.numeric(table(cellchat.plx.00d@idents))
# netVisual_circle(cellchat.plx.00d@net$count, vertex.weight = groupSize.plx.00d, weight.scale = T, 
#                  label.edge= F, title.name = "Number of interactions - 00d")
# netVisual_circle(cellchat.plx.00d@net$weight, vertex.weight = groupSize.plx.00d, weight.scale = T, 
#                  label.edge= F, title.name = "Interaction weights/strength - 00d")
# 
# # control 07d
# cellchat.plx.07d <- use_CellChat(counts.plx.07d, meta.plx.07d)
# cellchat.plx.07d <- aggregateNet(cellchat.plx.07d)
# groupSize.plx.07d <- as.numeric(table(cellchat.plx.07d@idents))
# netVisual_circle(cellchat.plx.07d@net$count, vertex.weight = groupSize.plx.07d, weight.scale = T, 
#                  label.edge= F, title.name = "Number of interactions - 07d")
# netVisual_circle(cellchat.plx.07d@net$weight, vertex.weight = groupSize.plx.07d, weight.scale = T, 
#                  label.edge= F, title.name = "Interaction weights/strength - 07d")
# 
# # control 28d
# cellchat.plx.28d <- use_CellChat(counts.plx.28d, meta.plx.28d)
# cellchat.plx.28d <- aggregateNet(cellchat.plx.28d)
# groupSize.plx.28d <- as.numeric(table(cellchat.plx.28d@idents))
# netVisual_circle(cellchat.plx.28d@net$count, vertex.weight = groupSize.plx.28d, weight.scale = T, 
#                  label.edge= F, title.name = "Number of interactions - 28d")
# netVisual_circle(cellchat.plx.28d@net$weight, vertex.weight = groupSize.plx.28d, weight.scale = T, 
#                  label.edge= F, title.name = "Interaction weights/strength - 28d")
# 
# dev.off()
# 
# # get the information of pathways
# p1 <- cellchat.ctrl.00d@netP$pathways # show all relevant pathways
# p2 <- cellchat.ctrl.07d@netP$pathways # show all relevant pathways
# p3 <- cellchat.ctrl.28d@netP$pathways # show all relevant pathways
# p4 <- cellchat.plx.00d@netP$pathways # show all relevant pathways
# p5 <- cellchat.plx.07d@netP$pathways # show all relevant pathways
# p6 <- cellchat.plx.28d@netP$pathways # show all relevant pathways
# Reduce(intersect, list(p1, p2, p3, p4, p5, p6)) # get the intersected pathways
# pathways.show <- c("JAM") # use this pathway as an example
# 
# # generate circle plots
# pdf("circle_plots.pdf")
# par(mfrow = c(2, 3), xpd=TRUE)
# #vertex.receiver <- c(1, 3) # select the target cell types
# netVisual_aggregate(cellchat.ctrl.00d, signaling = pathways.show, layout = "circle")
# netVisual_aggregate(cellchat.ctrl.07d, signaling = pathways.show, layout = "circle")
# netVisual_aggregate(cellchat.ctrl.28d, signaling = pathways.show, layout = "circle")
# netVisual_aggregate(cellchat.plx.00d, signaling = pathways.show, layout = "circle")
# netVisual_aggregate(cellchat.plx.07d, signaling = pathways.show, layout = "circle")
# netVisual_aggregate(cellchat.plx.28d, signaling = pathways.show, layout = "circle")
# dev.off()
# 
# 
# # generate the chord plots
# pdf("chord_plots.pdf")
# par(mfrow = c(2, 3), xpd=TRUE)
# #vertex.receiver <- c(1, 3) # select the target cell types
# netVisual_aggregate(cellchat.ctrl.00d, signaling = pathways.show, layout = "chord")
# netVisual_aggregate(cellchat.ctrl.07d, signaling = pathways.show, layout = "chord")
# netVisual_aggregate(cellchat.ctrl.28d, signaling = pathways.show, layout = "chord")
# netVisual_aggregate(cellchat.plx.00d, signaling = pathways.show, layout = "chord")
# netVisual_aggregate(cellchat.plx.07d, signaling = pathways.show, layout = "chord")
# netVisual_aggregate(cellchat.plx.28d, signaling = pathways.show, layout = "chord")
# dev.off()
# 
# 
# ## generate the heatmaps
# #pdf("heatmaps.pdf")
# #par(mfrow = c(2, 3), xpd=TRUE)
# ##vertex.receiver <- c(1, 3) # select the target cell types
# #netVisual_heatmap(cellchat.ctrl.00d, signaling = pathways.show, color.heatmap  = "Reds", 
# #                  title.name = paste0("Heatmap for ", pathways.show[1], " pathway - control 00d"))
# ##netVisual_heatmap(cellchat.ctrl.07d, signaling = pathways.show, color.heatmap  = "Reds", 
# ##                  title.name = paste0("Heatmap for ", pathways.show[1], " pathway - control 07d"))
# #netVisual_heatmap(cellchat.ctrl.28d, signaling = pathways.show, color.heatmap  = "Reds", 
# #                  title.name = paste0("Heatmap for ", pathways.show[1], " pathway - control 28d"))
# #netVisual_heatmap(cellchat.plx.00d, signaling = pathways.show, color.heatmap  = "Reds", 
# #                  title.name = paste0("Heatmap for ", pathways.show[1], " pathway - PLX5622 00d"))
# #netVisual_heatmap(cellchat.plx.07d, signaling = pathways.show, color.heatmap  = "Reds", 
# #                  title.name = paste0("Heatmap for ", pathways.show[1], " pathway - PLX5622 07d"))
# #netVisual_heatmap(cellchat.plx.28d, signaling = pathways.show, color.heatmap  = "Reds", 
# #                  title.name = paste0("Heatmap for ", pathways.show[1], " pathway - PLX5622 28d"))
# #dev.off()
# 
# 
# ## calculate the contribution of each ligand-receptor pair to the overall pathway
# #pdf("contributions.pdf")
# #par(mfrow = c(2, 3), xpd=TRUE)
# #netAnalysis_contribution(cellchat.ctrl.00d, signaling = pathways.show, 
# #                         title = paste0("Ligand-receptor contributions to ", pathways.show[1], " pathway - control 00d"))
# #netAnalysis_contribution(cellchat.ctrl.07d, signaling = pathways.show, 
# #                         title = paste0("Ligand-receptor contributions to ", pathways.show[1], " pathway - control 07d"))
# #netAnalysis_contribution(cellchat.ctrl.28d, signaling = pathways.show, 
# #                         title = paste0("Ligand-receptor contributions to ", pathways.show[1], " pathway - control 28d"))
# #netAnalysis_contribution(cellchat.plx.00d, signaling = pathways.show, 
# #                         title = paste0("Ligand-receptor contributions to ", pathways.show[1], " pathway - PLX5622 00d"))
# #netAnalysis_contribution(cellchat.plx.07d, signaling = pathways.show, 
# #                         title = paste0("Ligand-receptor contributions to ", pathways.show[1], " pathway - PLX5622 07d"))
# #netAnalysis_contribution(cellchat.plx.28d, signaling = pathways.show, 
# #                         title = paste0("Ligand-receptor contributions to ", pathways.show[1], " pathway - PLX5622 28d"))
# #dev.off()
# 
# ## bubble plot
# #pdf("bubble_plots.pdf")
# #par(mfrow = c(2, 3), xpd=TRUE)
# #netVisual_bubble(cellchat.ctrl.00d, sources.use = 1, targets.use = c(2, 3), remove.isolate = FALSE)
# #netVisual_bubble(cellchat.ctrl.07d, sources.use = 1, targets.use = c(2, 3), remove.isolate = FALSE)
# #netVisual_bubble(cellchat.ctrl.28d, sources.use = 1, targets.use = c(2, 3), remove.isolate = FALSE)
# #netVisual_bubble(cellchat.plx.00d, sources.use = 1, targets.use = c(2, 3), remove.isolate = FALSE)
# #netVisual_bubble(cellchat.plx.07d, sources.use = 1, targets.use = c(2, 3), remove.isolate = FALSE)
# #netVisual_bubble(cellchat.plx.28d, sources.use = 1, targets.use = c(2, 3), remove.isolate = FALSE)
# #dev.off()
# 
# # display the L-R pairs in chord plots
# pdf("chord_gene_diagrams.pdf")
# par(mfrow = c(2, 3), xpd=TRUE)
# netVisual_chord_gene(cellchat.ctrl.00d, sources.use = 1, targets.use = c(2, 3), lab.cex = 0.5,legend.pos.y = 30, 
#                      title.name = "Chord diagram - control 00d")
# #netVisual_chord_gene(cellchat.ctrl.07d, sources.use = 1, targets.use = c(2, 3), lab.cex = 0.5,legend.pos.y = 30, 
# #                     title.name = "Chord diagram - control 07d")
# netVisual_chord_gene(cellchat.ctrl.28d, sources.use = 1, targets.use = c(2, 3), lab.cex = 0.5,legend.pos.y = 30, 
#                      title.name = "Chord diagram - control 28d")
# netVisual_chord_gene(cellchat.plx.00d, sources.use = 1, targets.use = c(2, 3), lab.cex = 0.5,legend.pos.y = 30, 
#                      title.name = "Chord diagram - PLX5622 00d")
# netVisual_chord_gene(cellchat.plx.07d, sources.use = 1, targets.use = c(2, 3), lab.cex = 0.5,legend.pos.y = 30, 
#                      title.name = "Chord diagram - PLX5622 07d")
# netVisual_chord_gene(cellchat.plx.28d, sources.use = 1, targets.use = c(2, 3), lab.cex = 0.5,legend.pos.y = 30, 
#                      title.name = "Chord diagram - PLX5622 28d")
# dev.off()
