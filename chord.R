# Using this code, we can generate chord diagrams for all cell types in a dataset

# load libraries
library(CellChat)
library(patchwork)
library(Hmisc)
library(hash)
library(tidyverse)
library(circlize)
library(ggsci)
library(igraph)
library(gtools)
library(ComplexHeatmap)

# set the working directory
data.dir <- './'
setwd(data.dir)

# cell types
ct_list <- c('Macrophage', 'Astrocyte', 'Microglia')


# load data
{
  counts.ctrl.00d <- read.table("counts_control_00d.txt", header = T)
  counts.ctrl.00d <- normalizeData(counts.ctrl.00d)
  meta.ctrl.00d <- read.table("meta_control_00d.txt", header = T)
  
  counts.ctrl.07d <- read.table("counts_control_07d.txt", header = T)
  counts.ctrl.07d <- normalizeData(counts.ctrl.07d)
  meta.ctrl.07d <- read.table("meta_control_07d.txt", header = T)
  
  counts.ctrl.28d <- read.table("counts_control_28d.txt", header = T)
  counts.ctrl.28d <- normalizeData(counts.ctrl.28d)
  meta.ctrl.28d <- read.table("meta_control_28d.txt", header = T)
  
  counts.plx.00d <- read.table("counts_PLX5622_00d.txt", header = T)
  counts.plx.00d <- normalizeData(counts.plx.00d)
  meta.plx.00d <- read.table("meta_PLX5622_00d.txt", header = T)
  
  counts.plx.07d <- read.table("counts_PLX5622_07d.txt", header = T)
  counts.plx.07d <- normalizeData(counts.plx.07d)
  meta.plx.07d <- read.table("meta_PLX5622_07d.txt", header = T)
  
  counts.plx.28d <- read.table("counts_PLX5622_28d.txt", header = T)
  counts.plx.28d <- normalizeData(counts.plx.28d)
  meta.plx.28d <- read.table("meta_PLX5622_28d.txt", header = T)
}

# the function to use CellChat
use_CellChat <- function(counts, meta) {
  
  # counts : the gene expression matrix
  # meta   : the meta data matched with the gene expression matrix
  
  counts <- as.matrix(counts) # transform data frame into matrix
  cellchat <- createCellChat(object = counts, meta = meta, group.by = "Cell_type") # create a CellChat object
  levels(cellchat@idents) # check the cell types
  CellChatDB <- CellChatDB.mouse # set the database according to the organism
  dplyr::glimpse(CellChatDB$interaction) # show the structure of the database
  unique(CellChatDB$interaction$annotation) # check the annotation of interactions
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
    gene_pair1 <- df.net$ligand[x] # genes corresponding to the ligand
    genes1 <- capitalize(tolower(strsplit(gene_pair1, split = '_')[[1]]))
    gene_pair2 <- df.net$receptor[x] # genes corresponding to the receptor
    genes2 <- capitalize(tolower(strsplit(gene_pair2, split = '_')[[1]]))
    if (length(genes2) > 1 & genes2[2] == 'R2') {
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
}

id <- 0
df_list <- c() # the list to hold the CellChat object

# run CellChat for all datasets
{
  # control 00d
  net.ctrl.00d <- use_CellChat(counts.ctrl.00d, meta.ctrl.00d)
  id <- id + 1
  df_list[[id]] <- 'net.ctrl.00d'
  
  # control 07d
  net.ctrl.07d <- use_CellChat(counts.ctrl.07d, meta.ctrl.07d)
  id <- id + 1
  df_list[[id]] <- 'net.ctrl.07d'
  
  # control 28d
  net.ctrl.28d <- use_CellChat(counts.ctrl.28d, meta.ctrl.28d)
  id <- id + 1
  df_list[[id]] <- 'net.ctrl.28d'
  
  # PLX 00d
  net.plx.00d <- use_CellChat(counts.plx.00d, meta.plx.00d)
  id <- id + 1
  df_list[[id]] <- 'net.plx.00d'
  
  # PLX 07d
  net.plx.07d <- use_CellChat(counts.plx.07d, meta.plx.07d)
  id <- id + 1
  df_list[[id]] <- 'net.plx.07d'
  
  # PLX 28d
  net.plx.28d <- use_CellChat(counts.plx.28d, meta.plx.28d)
  id <- id + 1
  df_list[[id]] <- 'net.plx.28d'
}

# save the adjacency matrix and L-R categories
save_data <- function(net.df, ct_list, prefix) {
  
  # net.df  : the data frame to save the CCIs
  # ct_list : the list of all cell types
  # prefix  : the prefix of the files to output
  
  node_v1 <- vector() # the first node
  node_v2 <- vector() # the second node
  weight_v <- vector() # the edge weight
  
  for (i in 1:nrow(net.df)) {
    cci <- strsplit(net.df$interaction_name_2[i], split = ' ')[[1]] # split the interaction name
    ct.pair <- strsplit(net.df$cell_type_pair[i], split = '_')[[1]] # get the cell type pair
    ct1 <- ct.pair[1]
    ct2 <- ct.pair[2]
    prefix1 <- substr(ct1, start = 1, stop = 3)
    prefix2 <- substr(ct2, start = 1, stop = 3)
    node1 <- paste0(prefix1, '_', cci[1])
    id2 <- cci[length(cci)]
    weight <- net.df$prob[i]
    if (length(grep('\\+', id2)) > 0) {
      array <- strsplit(id2, split = '\\+')[[1]]
      node2 <- paste0(prefix2, '_', substr(array[1], 2, nchar(array[1])))
      node3 <- paste0(prefix2, '_', substr(array[2], 1, nchar(array[2]) - 1))
      node_v1 <- c(node_v1, node1, node1)
      node_v2 <- c(node_v2, node2, node3)
      weight_v <- c(weight_v, weight, weight)
    } else {
      node2 <- paste0(prefix2, '_', id2)
      node_v1 <- c(node_v1, node1)
      node_v2 <- c(node_v2, node2)
      weight_v <- c(weight_v, weight)
    }
  }
  
  if (length(node_v1) <= 0) {
    return()
  }
  
  g_df <- data.frame(node_v1, node_v2, weight_v)
  g <- graph.data.frame(g_df, directed = T)
  E(g)$weight <- g_df[[3]]
  adj <- get.adjacency(g, attr = 'weight')
  write.csv(x = as.data.frame(as.matrix(adj)), file = paste0(prefix, '_mat.csv'))
  
  Genes <- names(V(g))
  arrays <- strsplit(Genes, split = '_')
  ID <- sapply(arrays, function(x) {
    return(x[1])
  })
  category <- data.frame(Genes = Genes, ID = ID)
  write.csv(x = category, file = paste0(prefix, '_cat.csv'))
}

# convert the CCIs into the format to plot chord diagrams
for (i in 1:length(df_list)) {
  condition <- df_list[[i]] # the experiment condition
  prefix <- paste0(gsub(pattern = '^net.', replacement = '', df_list[[i]]))
  save_data(get(df_list[i][[1]]), ct_list, prefix)
}

# plot a chord diagram for a single datasete
genChord <- function(mat, cat, prefix) {
  
  # mat    : the adjacency matrix for the CCIs
  # cat   : the category for each gene, i.e., cell type
  # prefix : the prefix of the figure and legend files
  
  library(tidyverse)
  library(circlize)
  library(ggsci)
  library(igraph)
  library(gtools)
  library(ComplexHeatmap)
  
  # read data
  graph_adj <- read.csv(mat, row.names = 1)
  graph_module <- read.csv(cat, row.names = 1)
  g <- graph.adjacency(as.matrix(graph_adj), weighted = T)
  
  # set up color palette
  module_color <- pal_locuszoom()(7)
  # show_col(pal_locuszoom("default")(7))
  
  # Blue: #357ebdff
  # Brown: #7E6148FF
  # Turquoise: #4DBBD5FF
  # Yellow: #EEA236FF
  
  graph_module <- graph_module %>%
    mutate(color = as_factor(ID)) %>% 
    mutate(color = fct_recode(color, 
                              "#357ebdff" = "Mac",
                              "#7E6148FF" = "Mic",
                              # "#4DBBD5FF" = "turquoise.module",
                              "#EEA236FF" = "Ast"))
  
  # filter graph by weights
  g <- graph.adjacency(as.matrix(graph_adj), weighted = T)
  
  raw_edges <-  as.data.frame(cbind(get.edgelist(g), E(g)$weight)) %>%
    mutate(
      V1 = gsub('\\.', '-', V1),
      V2 = gsub('\\.', '-', V2),
      V3 = as.numeric(V3),
      V4 = 1
    )
  edges <- raw_edges %>%
    arrange(V3)
  
  # Normalize the weight score to quantiles
  # quartiles_weight <- quantcut(edges$V3, 4)
  # levels(quartiles_weight) <- c(1:4)
  # edges$V3 <- as.integer(quartiles_weight)
  
  
  nodes <-  unique(c(edges$V1, edges$V2))
  
  sectors <- sort(unique(c(raw_edges$V1, raw_edges$V2)))
  
  # diagram colors
  grid_col <- graph_module %>%
    dplyr::filter(Genes %in% nodes) %>% 
    dplyr::select(Genes, color) %>%
    mutate(color = as.character(color)) %>%
    deframe()
  
  #https://vuetifyjs.com/en/styles/colors/#material-colors
  col_fun = colorRamp2(range(edges$V3), c("#FFFDE7", "#013220"))
  
  # plot diagram
  png(
    paste(paste0(prefix, "_chord.png"), sep = ""),
    width = 3500,
    height = 3500,
    res = 300
  )
  
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(-0.15,0.2))
  circos.initialize(sectors, xlim = c(0, 1))
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.border = NA)
  
  # we go back to the first track and customize sector labels
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name = get.cell.meta.data("sector.index")
      xlim = get.cell.meta.data("xlim")
      this_node_text_color <- graph_module %>%
        dplyr::filter(Genes == sector.name) %>%
        pull(color) %>%
        as.character()
      
      circos.rect(
        xlim[1],
        0,
        xlim[2],
        1,
        col = this_node_text_color,
        border = NA
      )
      
      circos.text(
        mean(xlim),
        2,
        CELL_META$sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        col= this_node_text_color
      )
    },
    bg.border = NA
  ) 
  
  for (i in seq_len(nrow(edges))) {
    link <- edges[i,]
    circos.link(link[[1]],
                c(0, 1),
                link[[2]],
                c(0, 1),
                col = col_fun(link[[3]]),
                border = NA)
  }
  
  dev.off()
  
  # plot legend
  lgd <- Legend(title ="Score", col_fun = col_fun)
  grid.draw(lgd)
  
  png(
    paste(paste0(prefix, "_legend.png"), sep = ""),
    width = 1000,
    height = 1000,
    res = 300
  )
  grid.draw(lgd)
  dev.off()
  
}

# plot diagrams in batch
mat.list <- Sys.glob("*_mat.csv") # the list of adjacency matrix files
cat.list <- Sys.glob("*_cat.csv") # the list of the gene category files, i.e., cell types

if (length(mat.list != cat.list)) {
  stop("The number of adjacency matrices is not equal to that of gene category files.\n")
}

for (i in 1:length(mat.list)) {
  genChord(mat = mat.list[i], cat = cat.list[i], 
           prefix = strsplit(x = mat.list[i], split = '_')[[1]][1])
}