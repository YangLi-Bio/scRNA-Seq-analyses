# load libraries
library(clusterProfiler)
library(org.Mm.eg.db)

# read files into program
microglia.raw.data <- read.csv("Microglia_raw_data.csv")
rownames(microglia.raw.data) <- microglia.raw.data[,1]
microglia.raw.data <- microglia.raw.data[,-1]

macrophage.raw.data <- read.csv("Macrophage_raw_data.csv")
rownames(macrophage.raw.data) <- macrophage.raw.data[,1]
macrophage.raw.data <- macrophage.raw.data[,-1]

astrocyte.raw.data <- read.csv("Astrocyte_raw_data.csv")
rownames(astrocyte.raw.data) <- astrocyte.raw.data[,1]
astrocyte.raw.data <- astrocyte.raw.data[,-1]

# the function to transform gene symbols to ensembl ids
symbolToEnsembl <- function(m)
# m : the input matrix containing gene symbols
{
  symbols <- rownames(m) # get the gene names
  ensembls <- bitr(symbols, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
  o <- m[ensembls[,1],] # only retain the genes having correspondence to ensembl ids
  rownames(o) <- ensembls[,2] # use ensembl ids as row names
  return(o) # return the matrix with gene ensembl ids
}

# the function to generate metadata
matToMeta <- function(m, ct)
# m  : the expression matrix
# ct : the cell type
{
  cells <- colnames(m)
  cts <- rep(ct, length(cells)) # get the cell type
  cell.list <- strsplit(cells, split = "\\.|\\_")
  times <- sapply(cell.list, function(x) {
    t <- paste0(x[1], x[2])
    t <- gsub(t, pattern = "^X", replacement = "")
    return(t)
  }) # get the time
  treatments <- sapply(cell.list, function(x) x[3]) # get the treatments
  
  out <- data.frame(cells, cts, times, treatments) # build a list
  colnames(out) <- c("Cell", "Cell_type", "Time", "Treatment")
  
  return(out)
}

# convert gene symbols to ensembl ids
# microglia.ensembls <- symbolToEnsembl(microglia.raw.data)
# macrophage.ensembls <- symbolToEnsembl(macrophage.raw.data)
# astrocyte.ensembls <- symbolToEnsembl(astrocyte.raw.data)

# CellChat does not require using Ensembl ID
microglia.ensembls <- microglia.raw.data
macrophage.ensembls <- macrophage.raw.data
astrocyte.ensembls <- astrocyte.raw.data

# check whether the three cell types have the same number of genes expressed
#length(intersect(rownames(microglia.ensembls), rownames(macrophage.ensembls)))
#length(intersect(rownames(microglia.ensembls), rownames(astrocyte.ensembls)))
#length(rownames(microglia.ensembls))
a <- rownames(microglia.ensembls) == rownames(macrophage.ensembls)
a[which(a == F)]
b <- rownames(microglia.ensembls) == rownames(astrocyte.ensembls)
b[which(b == F)]

# re-arrange the order of genes
macrophage.ensembls <- macrophage.ensembls[rownames(microglia.ensembls),]
astrocyte.ensembls <- astrocyte.ensembls[rownames(astrocyte.ensembls),]

# get metadata
microglia.meta <- matToMeta(microglia.ensembls, "Microglia")
write.csv(microglia.meta, "Microglia_meta.csv")

macrophage.meta <- matToMeta(macrophage.ensembls, "Macrophage")
write.csv(macrophage.meta, "Macrophage_meta.csv")

astrocyte.meta <- matToMeta(astrocyte.ensembls, "Astrocyte")
write.csv(astrocyte.meta, "Astrocyte_meta.csv")

# merge the three datasets in a column-by-column manner
counts <- cbind(microglia.ensembls, macrophage.ensembls, astrocyte.ensembls)
#counts <- cbind(Gene = rownames(counts), counts)
meta <- rbind(microglia.meta, macrophage.meta, astrocyte.meta)

# save the merged dataset
#write.table(counts, "counts.txt", row.names = F, sep = "\t", quote = F)
#write.table(meta, "meta.txt", sep = "\t", quote = F)

# partition the data according to the time

meta.ctrl.00d <- meta[which(meta$Time == '00d' & meta$Treatment == 'Control'),]
meta.ctrl.07d <- meta[which(meta$Time == '07d' & meta$Treatment == 'Control'),]
meta.ctrl.28d <- meta[which(meta$Time == '28d' & meta$Treatment == 'Control'),]

counts.ctrl.00d <- counts[, meta.ctrl.00d$Cell]
counts.ctrl.07d <- counts[, meta.ctrl.07d$Cell]
counts.ctrl.28d <- counts[, meta.ctrl.28d$Cell]

#counts.ctrl.00d <- cbind(Gene = rownames(counts.ctrl.00d), counts.ctrl.00d)
#counts.ctrl.07d <- cbind(Gene = rownames(counts.ctrl.07d), counts.ctrl.07d)
#counts.ctrl.28d <- cbind(Gene = rownames(counts.ctrl.28d), counts.ctrl.28d)

write.table(counts.ctrl.00d, "counts_control_00d.txt", sep = "\t", quote = F)
write.table(counts.ctrl.07d, "counts_control_07d.txt", sep = "\t", quote = F)
write.table(counts.ctrl.28d, "counts_control_28d.txt", sep = "\t", quote = F)
write.table(meta.ctrl.00d, "meta_control_00d.txt", sep = "\t", quote = F, row.names = F)
write.table(meta.ctrl.07d, "meta_control_07d.txt", sep = "\t", quote = F, row.names = F)
write.table(meta.ctrl.28d, "meta_control_28d.txt", sep = "\t", quote = F, row.names = F)

meta.plx.00d <- meta[which(meta$Time == '00d' & meta$Treatment == 'PLX5622'),]
meta.plx.07d <- meta[which(meta$Time == '07d' & meta$Treatment == 'PLX5622'),]
meta.plx.28d <- meta[which(meta$Time == '28d' & meta$Treatment == 'PLX5622'),]

counts.plx.00d <- counts[, meta.plx.00d$Cell]
counts.plx.07d <- counts[, meta.plx.07d$Cell]
counts.plx.28d <- counts[, meta.plx.28d$Cell]

#counts.plx.00d <- cbind(Gene = rownames(counts.plx.00d), counts.plx.00d)
#counts.plx.07d <- cbind(Gene = rownames(counts.plx.07d), counts.plx.07d)
#counts.plx.28d <- cbind(Gene = rownames(counts.plx.28d), counts.plx.28d)

write.table(counts.plx.00d, "counts_PLX5622_00d.txt", sep = "\t", quote = F)
write.table(counts.plx.07d, "counts_PLX5622_07d.txt", sep = "\t", quote = F)
write.table(counts.plx.28d, "counts_PLX5622_28d.txt", sep = "\t", quote = F)
write.table(meta.plx.00d, "meta_PLX5622_00d.txt", sep = "\t", quote = F, row.names = F)
write.table(meta.plx.07d, "meta_PLX5622_07d.txt", sep = "\t", quote = F, row.names = F)
write.table(meta.plx.28d, "meta_PLX5622_28d.txt", sep = "\t", quote = F, row.names = F)