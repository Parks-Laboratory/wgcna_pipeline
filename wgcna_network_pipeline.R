#------------------------------------------------------------------------
# Title: Run WGCNA Wrapper Functions
# Date: May 23 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS: Functions for
  # Process SQL Data
  # Choosing Soft Threshold Power: Pick Soft Power
  # WGCNA wrapper: Run WGCNA 
  # Format for Cytoscape: Format Cytoscape
  # Clinical Traits Correlation: Clinical Trait Cor


# NOTES:
#------------------------------------------------------------------------
library(plyr)
library(stringr)
library(dplyr)
library(WGCNA)
try_default( source("E:/MICROARRAY DATA/wgcna_pipeline/wgcna_functions.R"), default = NULL, quiet = TRUE )


# PRELIMINARY
#-------------------------------------------------------------------------------------------------------------------


# Process SQL Data
#------------------------------------------------------------------------

# Parameters
  # sql_data: named list of datasets that came from SQL query calls 
  # use_col: (string) name of the column representing the expression values
  # output_folder: output folder

# Notes:
  # sql_data should have 3 columns: 
    # strain, which will go long
    # probeid, which will go wide
    # colname saved into the variable use_col, which will go in value cells

process_sql_data <- function(sql_data, use_col, outlier_cutoff = 100, output_folder = ""){
  
  # Reformat Data
  #------------------------------------------------------------------------
  
  # format data for networks: strains long probes wide
  expr_data <- lapply(sql_data, function(d){
    data <- as.data.frame(dcast(d, strain ~ probe_id, value.var = use_col))
    rownames(data) <- data$strain
    data <- dplyr::select(data, -strain)
    colnames(data) <- paste0("P", colnames(data))
    return(data)
  })
  
  
  # QC: Missing Data and Outliers
  #------------------------------------------------------------------------
  
  # apply data checks
  expr_data_after <- lapply(1:length(expr_data), function(i){
    
    datExpr0 <- expr_data[[i]]
    
    # check for excessive missing values and identify outlier microarrays
    cat("Checking for missing values")
    
    # first check for genes and samples with too many missing values
    gsg <- goodSamplesGenes(datExpr0, verbose = 3)
    
    # remove offending genes & samples from the data set
    if(!gsg$allOK) {
      
      if(sum(!gsg$goodGenes) > 0) printFlush(paste('Removing genes:', paste(names(datExpr0)[!gsg$goodGenes], collapse = ', ')))
      if(sum(!gsg$goodSamples) > 0) printFlush(paste('Removing samples:', paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ', ')))
      
      # Remove offending genes and samples from the data
      datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
      
    }
    
    # check for outlier samples
    cat( paste('Checking for outliers: cutoff =', outlier_cutoff, '\n') )
    
    # cluster the samples to detect outliers: set cutoff point at 100
    sampleTree <- flashClust::flashClust(dist(datExpr0), method = 'average')
    cut <- outlier_cutoff
    
    par(cex = 0.6, mar = c(0, 4, 2, 0))
    pdf(paste0(output_folder, "/sample_cluster_for_outliers_", names(expr_data)[i], ".pdf"), height = 9, width = 12)
    plot(sampleTree, main = 'Sample clustering to detect outliers', xlab='', cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    abline(h = cut, col = "red")
    dev.off()
    
    # find the clusters underneath the line and keep specified samples
    clust <- cutreeStatic(sampleTree, cutHeight = cut, minSize = 10 )
    keep <- names(which.max(table(clust)))
    datExpr <- datExpr0[clust == keep, ]
    
    # return results
    return(datExpr)
  })
  names(expr_data_after) <- names(expr_data)
  
  # QC: express data dimensions prior to fix
  message("Expression data dimensions before fix")
  print( lapply(expr_data, dim) )
  
  # QC: express data dimensions prior to fix
  message("Expression data dimensions after fix")
  print( lapply(expr_data_after, dim) )
  
  # change name
  expr_data <- expr_data_after
  
  # extract strain names
  strain_names <- lapply(expr_data, rownames)
  
  
  # Save Results into Rdata File
  #------------------------------------------------------------------------
  
  # save the results into an expression file
  message( paste("Saving Expression Data in", output_folder, "as network_expr_data.Rdata") )
  save(expr_data, strain_names, file = file.path(output_folder, "network_expr_data.Rdata"))
  
}


# Choose the soft-thresholding power
#------------------------------------------------------------------------

# Parameters
  # datExpr: expression data
  # output_folder: folder for outputs
  # output_suffix: output file suffix

pick_soft_power <- function(datExpr, output_folder = "", output_suffix = ""){
  
  # choose a set of powers
  powers <- c( 1:10, seq(from = 12, to = 20, by = 2) )
  
  # call the network topology analysis function: what should we set soft power to?
  soft_threshold <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  
  # plot the results
  par(mfrow = c(1,2))
  cex1 <- 0.9
  pdf(file.path(output_folder, paste0("scale_independence_and_mean_connectivity", output_suffix, ".pdf")), width=9, height = 5 )
  plot(soft_threshold$fitIndices[,1], -sign(soft_threshold$fitIndices[,3])*soft_threshold$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
  text(soft_threshold$fitIndices[,1], -sign(soft_threshold$fitIndices[,3])*soft_threshold$fitIndices[,2], labels = powers, cex = cex1, col = "red")
  # this line corresponds to using a R^2 cutoff of h
  abline(h = 0.80, col = "red")
  
  # mean connectivity as a function of the soft-thresholding power
  plot(soft_threshold$fitIndices[,1], soft_threshold$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
  text(soft_threshold$fitIndices[,1], soft_threshold$fitIndices[,5], labels = powers, cex = cex1, col = "red")
  dev.off()
  
}



# PIPELINE
#-------------------------------------------------------------------------------------------------------------------


# Run WCGNA
#------------------------------------------------------------------------

# Function to run wcgna
# Parameters:
  # expr_data: expression data 
  # probeid: probe ids of the expression data
  # strains: name of the strains (rows)
  # data_id: identifier for data - will name the output files accordingly
  # affy: data frame that contains a column named probeid that matches the probeids from the expr data, gene_symbol, entrez_gene used to run GOenrichments; if no entrez_gene data, make it a column of NAs
  # clinical_traits: data frame of clinical traits; may leave empty (best if there are few modules <20 otherwise do it outside of pipeline)
  # softpower: softpower parameter (run pick_soft_power to find out)
  # min_mod_size: minimum size of the modules
  # deep_split: parameter for trimming dendrogram
  # ME_cor_thres: correlation threshold, merge modules of their eigengenes have correlation > 1 - ME_cor_thres
  # output_folder: folder for outputs
  # output_genetree: output genes/modules in order to plot clusters - used for clusters
# Returns: nothing; saves output

run_wgcna <- function(expr_data, probeid, strains, data_id = "", affy, clinical_traits, 
                      soft_power, min_mod_size = 30, deep_split = 3, ME_cor_thres = 0.1, 
                      output_folder = "", output_genetree = FALSE){
  
  # cluster genes
  genes <- cluster_genes(datExpr = expr_data, soft_power = soft_power, min_module_size = min_mod_size, deep_split = deep_split, output_folder = output_folder, output_suffix = paste0("_", data_id))
  
  # cluster modules
  modules <- cluster_modules(datExpr = expr_data, genes$geneTree, genes$dynamicColors, soft_power = soft_power, ME_cor_thres = ME_cor_thres, output_folder = output_folder, output_suffix = paste0("_", data_id))
  
  # add strain names to the ME data
  MEs <- output_MEs(modules, strains, output_folder, paste0(data_id, "_"))
  
  # format networks into df
  networks <- network_to_df(network_modules = modules, probeid = probeid, affy = affy, output_folder = output_folder, output_prefix = paste0(data_id, "_"))
  
  # intiate list
  network_data <- list(id = data_id, TOM = genes$TOM, geneTree = genes$geneTree, MEs = MEs, networkDF = networks, GOenrich = "GO enrichments not executed")
  
  # run GOenrichments
  if( any(str_detect(colnames(networks), "entrez")) ){
    
    # run GOenrichment
    GOenr <- runGO(networkdf = subset(networks, !is.na(entrez_id)), output_folder = output_folder, output_prefix = paste0(data_id, "_"))
    
    # organize output and return
    network_data$GOenrich <- GOenr

  } 

  # save output
  message( paste0("Saving Network Data in ", output_folder, " as ", data_id, "_networkConstruction.Rdata") )
  save(network_data, file = file.path(output_folder, paste0(data_id, "_networkConstruction.Rdata")))
  
  # run clinical trait information
  if( !missing(clinical_traits) ){
    
    # subset clinical traits
    sub_clinical_traits <- subset(clinical_traits, strain %in% MEs$strain)
    
    # run clinical trait correlations
    png(file.path(output_folder, paste0("clinical_trait_cor_", data_id, ".png")), height = 800, width = 1400)
    clinical_trait_cor(MEs, sub_clinical_traits)
    dev.off()

  }
  
  # output gene tree (if running on clusters)
  if(output_genetree){
    save(genes, modules, file = file.path(output_folder, paste0(data_id, "_gene_module_data.Rdata")))
  }
}


# Format for Cytoscape
#------------------------------------------------------------------------

# Function to format network data into formats for cytoscape
# Parameters: 
  # network_data: output from run_wcgna; may want to edit networkdf object to add additional attributes for cytoscape
  # expr_data: original expression data
  # modules: vector of module colors to include
  # threshold: threshold for TOM value for inclusion into the edges[0-1]
  # include_col: names of other columns in the network df to include (usually for coloring, etc)
  # output_folder: folder for outputs
  # output_prefix: output file prefix
# Returns: nothing

format_cytoscape <- function(network_data, expr_data, modules, threshold, include_col, output_folder = "", output_suffix = ""){
  
  # subset network data to specified module, removing duplicates
  sub_network <- network_data$networkDF %>% 
    subset(clusterL %in% modules) %>% 
    dplyr::select(probeid, clusterL) %>% 
    distinct()
  
  # combine probe id with gene symbols
  sub_network <- merge(sub_network, affy, "probeid", all.x = TRUE)
  
  # subset TOM data (fix if adding a prefix P to probe ids)
  prefix <- unique( str_sub(colnames(expr_data), 1, 1) )
  if( length(prefix) == 1 & prefix == "P"){
    keepTOM <- str_replace(colnames(expr_data), "^P", "") %in% sub_network$probeid 
    namesTOM <- str_replace(colnames(expr_data)[keepTOM], "^P", "")
  } else{
    keepTOM <- colnames(expr_data) %in% sub_network$probeid
    namesTOM <- colnames(expr_data)[keepTOM]
  }

  modTOM <- network_data$TOM[which(keepTOM), which(keepTOM)]
  dimnames(modTOM) <- list(namesTOM, namesTOM)
  
  # reorder network
  network <- merge(data.frame(probeid = colnames(modTOM)), sub_network, "probeid")
  
  # export the network into edge and node list files Cytoscape can read
  cyt <- exportNetworkToCytoscape(modTOM,
                                  edgeFile = file.path(output_folder, paste("CytoscapeInput-edges-", paste(modules, collapse=" - "), "-", output_suffix, ".txt", sep="")), 
                                  nodeFile = file.path(output_folder, paste("CytoscapeInput-nodes-", paste(modules, collapse=" - "), "-", output_suffix, ".txt", sep="")), 
                                  weighted = TRUE, 
                                  threshold = threshold,
                                  nodeNames = network$probeid, 
                                  altNodeNames = network$gene_symbol, 
                                  nodeAttr = dplyr::select(network, clusterL, one_of(include_col)))
  
}


# Clinical Traits Correlation
#------------------------------------------------------------------------

# Function to correlate module eigengenes with clinical traits
# Parameters: (rows of parameters must match)
  # mod_eigengenes: df of module eigengenes with strain column, sorted by strain
  # clinical_traits: df of clinical traits with strain column, sorted by strain
# Returns: heatmap of correlation table

clinical_trait_cor <- function(mod_eigengenes, clinical_traits){
  
  # sort data by strain info
  mod_eigengenes <- arrange(mod_eigengenes, strain)
  clinical_traits <- arrange(clinical_traits, strain)
  
  # throw error if the strains are not correct
  stopifnot( all(mod_eigengenes$strain == clinical_traits$strain) )

  # remove strain info
  MEs <- dplyr::select(mod_eigengenes, -strain)
  CTs <- dplyr::select(clinical_traits, -strain)
  
  # correlate eigengenes with trait data
  mod_correlate <- bicorAndPvalue(MEs, CTs)
  mod_bicor <- mod_correlate$bicor
  mod_p <- mod_correlate$p
  
  # generate visualization
  textMatrix <- paste0(signif(mod_bicor, 2), "\n(",signif(mod_p, 1), ")")
  dim(textMatrix) <- dim(mod_bicor)
  hm <- labeledHeatmap(Matrix = mod_bicor,
                 xLabels = colnames(CTs),
                 yLabels = colnames(MEs),
                 ySymbols = colnames(MEs) %>% str_replace("ME", ""),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = TRUE,
                 cex.text = 0.9,
                 xLabelsAngle = 30,
                 yColorWidth = 0.01,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  
  # return visualization
  return(hm)
  
}
  
  