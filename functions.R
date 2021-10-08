#########
# Load useful functions, do not print in the final report
#########


# invisibly calculate a value
quiet <- function(x) {
  sink(tempfile()) # Generate temporary files and directories
  on.exit(sink()) # sink : diverts R output to a connection, e.g., txt and csv, as parameter
  # on.exit : records the expression given as its argument as needing to be
  # executed when the current function exits
  invisible(force(x)) # force : forces the evaluation of a formal argument
  # Use invisible in place of return in a function when the returned output should not be printed
}


# point size function from test datasets
x <- c(0, 90, 124, 317, 1000, 2368, 3005, 4816, 8298, 50000, 500000, 5000000)
y <- c(1, 1, 0.89, 0.33, 0.30, 0.25, 0.235, 0.205, 0.18, 0.1, 0.1, 0.1)
get_point_size <- approxfun(x, y) # The function approxfun returns a function performing
# (linear or constant) interpolation of the given data points


# generate the UMAP plot for expression of a gene
Plot.GeneUMAP <- function(object = my.object, gene.name = NULL, pt_size = 0.5){
  tmp.gene.expression <- GetAssayData(object) # the expression matrix
  tmp.dim <- as.data.frame(object@reductions$umap@cell.embeddings) # two-column embedding axises
  tmp.MatchIndex <- match(colnames(tmp.gene.expression), rownames(tmp.dim))
  # match returns a vector of the positions of (first) matches of its first argument in its second
  tmp.dim <- tmp.dim[tmp.MatchIndex, ]
  tmp.gene.name <- paste0("^", gene.name, "$")
  tmp.One.gene.value <- tmp.gene.expression[grep(tmp.gene.name, rownames(tmp.gene.expression)), ]
  # the expression profile of one gene
  tmp.dim.df <- cbind.data.frame(tmp.dim, Gene = tmp.One.gene.value) # three column data frame
  g <- ggplot(tmp.dim.df, aes(x = UMAP_1, y= UMAP_2, color = Gene))
  g <- g + geom_point(stroke = pt_size, size = pt_size) + scale_color_gradient(low = "grey", high = "red")
  # Use the stroke aesthetic to modify the width of the border
  g <- g + theme_classic() + labs(color = paste0(gene.name, "\nexpression\nvalue")) + coord_fixed(1)
  # A classic-looking theme, with x and y axis lines and no gridlines
  # labs : Change axis labels and legend titles
  # coord_fixed : A fixed scale coordinate system forces a specified ratio between the physical
  #     representation of data units on the axes
  g
}


# similar to Plot.GeneUMAP
Plot.GeneTSNE <- function(object = my.object, gene.name = NULL, pt_size = 0.5) {
  tmp.gene.expression <- object@assays$SCT@data
  tmp.dim <-as.data.frame(object@reductions$tsne@cell.embeddings)
  tmp.MatchIndex <- match(colnames(tmp.gene.expression), rownames(tmp.dim))
  tmp.dim <- tmp.dim[tmp.MatchIndex, ]
  tmp.gene.name <- paste0("^", gene.name,"$")
  tmp.One.gene.value <- tmp.gene.expression[grep(tmp.gene.name,rownames(tmp.gene.expression)), ]
  tmp.dim.df <- cbind.data.frame(tmp.dim, Gene=tmp.One.gene.value)
  g <- ggplot(tmp.dim.df, aes(x = tSNE_1, y = tSNE_2, color = Gene))
  g <- g + geom_point(stroke = pt_size, size=pt_size) + scale_color_gradient(low = "grey", high = "red")
  g <- g + theme_classic() + labs(color = paste0(gene.name, "\nexpression\nvalue")) + coord_fixed(1)
  g
}


# generate 2D plot for cell clusters
Plot.cluster2D <- function(object = combined, reduction.method = "umap",
    customized = T, pt_size=1, txt = "Cell type", ...) {
  
  my.plot.all.source <- cbind.data.frame(Embeddings(object, reduction = reduction.method),
                                       Cell_type = Idents(object))
  
  tmp.celltype <- levels(unique(my.plot.all.source$Cell_type))
  p.cluster <- ggplot(my.plot.all.source,
                      aes(x = my.plot.all.source[, 1], y = my.plot.all.source[, 2])) +
                      xlab(colnames(my.plot.all.source)[1]) +
                      ylab(colnames(my.plot.all.source)[2])
  p.cluster <- p.cluster + geom_point(stroke = pt_size, size = pt_size,
                                aes(col = my.plot.all.source[, "Cell_type"]))
                                # Use the stroke aesthetic to modify the width of the # border
  
  p.cluster <- p.cluster + guides(colour = guide_legend(override.aes = list(size=5)))
  # guides can be used to modify the legend characteristics
  
  if (length(tmp.celltype) >= 5) {
    p.cluster <- p.cluster + scale_colour_manual(name = paste(txt, ":(Cells)", sep = ""), 
                                                 values = as.character(palette36.colors(36)[-2][1:length(tmp.celltype)]),
                                                 breaks = tmp.celltype,
                                                 labels = paste0(tmp.celltype, 
                                                 ":(",as.character(summary(my.plot.all.source$Cell_type)),")"))
  } else if (length(tmp.celltype) < 5) {
    p.cluster <- p.cluster + scale_colour_manual(name = paste(txt, ":(Cells)", sep = ""), 
                                                 values = brewer.pal(4, "Spectral")[c(2, 1, 3, 4)],
                                                 breaks = tmp.celltype,
                                                 labels = paste0(tmp.celltype, 
                                                 ":(",as.character(summary(my.plot.all.source$Cell_type)),")"))
  } else {
    p.cluster <- p.cluster + scale_colour_manual(name = paste(txt, ":(Cells)", sep = ""), 
                                                 values = brewer.pal(5,"Spectral")[c(1,5)],
                                                 breaks = tmp.celltype,
                                                 labels = paste0(tmp.celltype, 
                                                 ":(",as.character(summary(my.plot.all.source$Cell_type)),")"))
  }
  
  
  # + labs(col="cell type")           
  p.cluster <- p.cluster + theme_classic() 
  p.cluster <- p.cluster + coord_fixed(ratio = 1)
  p.cluster
}


##############################
# define a fucntion for reading in 10X hd5f data and cell gene matrix by input (TenX) or (CellGene)
# & and && : The shorter form performs elementwise comparisons in much the same way as arithmetic operators.

read_data <- function(x = NULL, read.method = NULL, sep="\t", ...) {
  if(!is.null(x)) {
    if(!is.null(read.method)) {
      if(read.method !="TenX.h5" && read.method != "CellGene" && read.method != "TenX.folder") {
        stop("wrong 'read.method' argument, please choose 'TenX.h5','TenX.folder', or 'CellGene'!")}
        
      if(read.method == "TenX.h5") {
        tmp_x <- Read10X_h5(x)
        return(tmp_x)
      } else if (read.method =="TenX.folder") {
        all_files <- list.files(getwd())
        barcode_file <- grep("barcodes", all_files)
        matrix_file <- grep("matrix", all_files)
        gene_file <- grep("genes", all_files)
        feature_file <- grep("features", all_files)
        
        #Check users upload single zipped file, by counting detected filename, if
        # less than 3 we think users uploads zipped file
        if((length(barcode_file) + length(matrix_file) + length(gene_file) +
            length(feature_file)) < 3) {
          dir.create("tmp", showWarnings = F)
          
          if (file_ext(x) == "7z") { # file_ext returns the file (name) extensions
          # (excluding the leading dot)
            try(system(paste("7za x", x, "-aoa -otmp")), silent = T)
            # try evaluates an expression and traps any errors that occur during the evaluation
          }
          
          try(system(paste("unzip -o", x, "-d tmp")), silent = T) # gip format
          try(system(paste("tar xzvf", x, "--directory tmp")), silent = T) # tar.gz format
          
          # check if the file is gz instead of tar.gz
          max_file <- which.max(file.info(list.files("tmp", full.names = T, recursive = T))[, 1])
          # choose the maximum file
          # file.info : extract information about files on the user's file systems including
          #     size, is.dir, mode, and gid, etc. as columns, and files as rows
          this_files <- list.files("tmp", full.names = T, recursive = T)[max_file] # the maximum file
          if(is.na(this_files) || file_ext(this_files) == "tar" || length(this_files) == 0) {
            system("rm -R tmp/*")
            this_filename <- gsub(".gz","", basename(x))
            # basename : removes all of the path up to and including the last path separator
            try(system(paste("gunzip -c ", x, " > tmp/", this_filename, sep = "")), silent = T)
            max_file <- which.max(file.info(list.files("tmp", full.names = T, recursive = T))[, 1])
            this_files <- list.files("tmp", full.names = T, recursive = T)[max_file]
            this_delim <- reader::get.delim(this_files)
            # The read.delim function is typically used to read in delimited text
            # returns character of the most likely delimiter
            tmp_z <- tryCatch(read.delim(paste0(this_files), header = T, row.names = NULL,
                            check.names = F, sep = this_delim), error = function(e) 0)
            # tryCatch : Evaluates an expression with the possibility to catch exceptions
            # read.delim : Reads a file in table format and creates a data frame
            #       from it, with cases corresponding to lines and variables to fields in the file.
            # check.names : ensure that they are syntactically valid variable names
            # error : Do this if an error is caught
            upload_type <<- "CellGene"
            return(tmp_z)
          }
          
          max_file <- which.max(file.info(list.files("tmp", full.names = T, recursive = T))[, 1])
          this_files <- list.files("tmp", full.names = T, recursive = T)[max_file]
          
          
          # in case folder contains 10X files
          tmp_x <- tryCatch(Read10X(gsub(basename(this_files), "", this_files)), error = function(e) 0)
          
          if (typeof(tmp_x) == "S4") {
            system("rm -R tmp/*")
            return(tmp_x)
          } else if (file_ext(this_files) == "h5" || file_ext(this_files) == "hdf5") {
            tmp_y <- tryCatch(Read10X_h5(this_files), error = function(e) 0)
            upload_type <<- "TenX.h5"
            system("rm -R tmp/*")
            
            return(tmp_y)
          } else {
            this_delim <- reader::get.delim(this_files)
            tmp_z <- tryCatch(read.delim(paste0(this_files), header = T, row.names = NULL,
                        check.names = F, sep = this_delim), error = function(e) 0)
            upload_type <<- "CellGene"
            system("rm -R tmp/*")
            return(tmp_z)
          }
          
        }
        
        # remove the keywords in file names
        tryCatch(file.rename(all_files[barcode_file], paste("barcodes", gsub(".*barcodes",
                                                            "", all_files[barcode_file]), sep = "")),
                                                            error = function(e) 0) # rename files
        tryCatch(file.rename(all_files[matrix_file], paste("matrix", gsub(".*matrix", "", 
                                                            all_files[matrix_file]), sep = "")),
                                                            error = function(e) 0)
        tryCatch(file.rename(all_files[gene_file], paste("genes", gsub(".*genes", "", all_files[gene_file]),
                                                            sep = "")), error = function(e) 0)
        tryCatch(file.rename(all_files[feature_file], paste("features", gsub(".*features", "",
                                                            all_files[features]), sep = "")),
                                                            error = function(e) 0)
        
        tmp_x <- tryCatch(Read10X(getwd()), error = function(e) {
          all_files <- list.files(getwd())
          barcode_file <- grep("barcodes", all_files)
          matrix_file <- grep("matrix", all_files)
          gene_file <- grep("genes", all_files)
          feature_file <- grep("features", all_files)
          try(system(paste("gunzip", (all_files[barcode_file]))), silent = T)
          try(system(paste("gunzip", (all_files[matrix_file]))), silent = T)
          try(system(paste("gunzip", (all_files[gene_file]))), silent = T)
          try(system(paste("gunzip", (all_files[feature_file]))), silent = T)
          try(system(paste("unzip", (all_files[barcode_file]))), silent = T)
          try(system(paste("unzip", (all_files[matrix_file]))), silent = T)
          try(system(paste("unzip", (all_files[gene_file]))), silent = T)
          try(system(paste("unzip", (all_files[feature_file]))), silent = T)
        })
        
        tmp_x <- tryCatch(Read10X(getwd()), error = function(e) {
          0
        })
        
        return(tmp_x)
      } else if (read.method == "CellGene"){ # read in cell * gene matrix, if there is error report, back to 18 line to run again.
        tmp_x < -read.delim(x, header = T, row.names = NULL, check.names = F, sep=sep,...)
        
        return(tmp_x)
      }
    }
  } else {stop("Missing 'x' argument, please input correct data")}
}


# convert the first letter to upper case
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  
  return (x)
}



## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot <- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin )
          
  return(p)
}


## extract the max value of the y axis
extract_max <- function(p) {
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  # ggplot_build() takes the plot object, and performs all steps necessary to produce
  # an object that can be rendered. This function outputs two pieces: a list of
  # data frames (one for each layer), and a panel object, which contain all
  # information about axis limits, breaks etc.
  
  return(ceiling(ymax))
}


## main function
StackedVlnPlot <- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  # unit : This function creates a unit object --- a vector of unit values. A unit value
  # is typically just a single numeric value with an associated unit.
  
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, feature = x, ...))
  # generate a violin plot for each feature
  # purrr::map() is a function for applying a function to each element of a list
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  
  return(p)
}
