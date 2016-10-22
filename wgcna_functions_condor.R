#------------------------------------------------------------------------
# Title: Run WGCNA Functions
# Date: May 23 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS: Functions for
  # Network Genes: Cluster Genes
  # Merge Modules: Cluster Modules
  # Edit Module Output to DF: Network to DF
  # Output Module Eigengenes: Output MEs
  # GSEA: Run GO


# NOTES:
#------------------------------------------------------------------------


# Run network (functions)
#------------------------------------------------------------------------

# Function that
  # 1) computes bicor matrix
  # 2) generates adjacency/TOM matrix to input into cluster algorithm
  # 3) hierarchical clustering using "average" linkage
  # 4) trim tree based on minModSize and deepSplit
  # 5) outputs plots
# Params:
  # datExpr: expression matrix
  # soft_power: power for correlations
  # min_module_size: minimum cluster size
  # deep_split: 0-4, rough control over sensitivity to cluster splitting; higher values = smaller/more clusters
  # output_folder: folder for outputs
  # output_suffix: output file suffix
# Returns: list of 
  # geneTree: dendogram of tree
  # dynamicColors: colors indicating the modules of the data

cluster_genes <- function(datExpr, soft_power, min_module_size, deep_split, output_folder = "", output_suffix = ""){
  
  # calculate biweight midcorrelation instead of default Pearson correlation in adjacency
  bicorMatrix <- bicor(datExpr, use = "pairwise.complete.obs")
  
  # generate adjacency matrix from bicorMatrix
  adjacency <- adjacency.fromSimilarity(bicorMatrix, type = "unsigned", power = soft_power)
  
  # transform adjacency into topological overlap matrix to minimize effects of noise and spurious associations
  dissTOM <- TOMdist(adjacency)
  
  # cluster genes using dissTOM
  geneTree <- flashClust::flashClust(as.dist(dissTOM), method = "average")
  
  # determine what parameter of deep split gives the best modules 
  # mColorh <- lapply(1:4, function(ds){
  #   tree <- cutreeHybrid(dendro = geneTree, pamStage = FALSE, minClusterSize = (min_module_size - 3*ds), cutHeight = 0.99, deepSplit = ds, distM = dissTOM)
  #   return( labels2colors(tree$labels) )
  # })
  # mColorh <- Reduce(cbind, mColorh)
  
  # print dynamic tree cut results
  # png(file.path(output_folder, paste0("deep_split_module_options", output_suffix, ".png")), height = 500, width = 700)
  # plotDendroAndColors(geneTree, mColorh, paste("deepSplit=", 1:4), main = "", dendroLabels = FALSE)
  # dev.off()
  
  # trim the tree based on input parameters
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = deep_split, pamRespectsDendro = FALSE, minClusterSize = min_module_size)
  
  # convert numeric labels into colors
  dynamicColors <- labels2colors(dynamicMods)
  
  # plot the final results
  # png(file.path(output_folder, paste0("gene_dendrogram_and_module_colors", output_suffix, ".png")), width = 700, height = 500)
  # plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE,  hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
  # dev.off()
  
  return(list(geneTree = geneTree, dynamicColors = dynamicColors, TOM = dissTOM))
}


# Function that
  # 1) computes the eigengenes and their dissimilarity
  # 2) generates dissimilarity using bicor
  # 3) hierarchical clusters using "average" linkage
  # 4) merge modules together based on threshold 
  # 5) outputs plots
# Params:
  # datExpr: expression matrix
  # tissue: tissue source of data (for output names)
  # geneTree: dendrogram of the gene clusters
  # soft_power: power for correlations
  # ME_cor_thres: threshold for merging modules, equals 1 - desired cor(modules)
  # output_folder: folder for outputs
  # output_suffix: output file suffix
# Returns: list of 
  # eigengenes & their variability, module groups, eigengenes of merged modules

cluster_modules <- function(datExpr, geneTree, dynamicColors, soft_power, ME_cor_thres, output_folder = "", output_suffix = ""){
  
  # calculate eigengenes
  MEList <- moduleEigengenes(datExpr, colors = dynamicColors, softPower = soft_power, excludeGrey = TRUE)
  MEvar <- MEList$varExplained
  MEs <- MEList$eigengenes
  
  # calculate dissimilarity of module eigengenes
  MEDiss <- 1 - bicor(MEs, use = "pairwise.complete.obs")
  
  # cluster modules whith similar expression profiles
  METree <- flashClust::flashClust(as.dist(MEDiss), method = "average")
  
  # plot module clusters
  # png(file.path(output_folder, paste0("clustering_module_eigengenes", output_suffix, ".png")), width = 700, height = 500)
  # plot(METree, main = "Clustering of module eigengenes", xlab = "")
  # abline(h = ME_cor_thres, col = "red")
  # dev.off()
  
  # merge modules together
  merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = ME_cor_thres, verbose = 3)
  mergedColors <- merge$colors
  
  # obtain eigengenes of merged modules
  mergedMEs <- dplyr::select(merge$newMEs, -matches("MEgrey"))
  
  # print merged modules
  # png(file.path(output_folder, paste0("clustering_after_merge", output_suffix, ".png")), height = 500, width = 700)
  # plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
  # dev.off()
  
  # return results
  return(list(MEs = MEs, MEvar = MEvar, mergedColors = mergedColors, mergedMEs = mergedMEs))
}


# Format network to a data frame
#------------------------------------------------------------------------

# Function that saves the network data in a csv
# Parameters:
  # network_modules: output from cluster_modules
  # probe_ids: probe ids of original expression data (in the specified order)
  # affy: data frame that contains a column named probeid that matches the probeids from the expr data, gene_symbol, entrez_gene (codes); if no entrez_gene data, make it a column of NAs
  # output_folder: folder for outputs
  # output_prefix: output file prefix
# Returns: data frame with columns
  # probeid
  # clusterL (cluster color)
  # gene_symbol
  # entrez_id

network_to_df <- function(network_modules, probeid, affy, output_folder = "", output_prefix = ""){
  
  prefix <- unique(str_sub(probeid, 1, 1))
  if( all(length(prefix) == 1 & prefix == "P") ) probe_id <- str_replace(probeid, "^P", "") else probe_id <- probeid
  
  # initalize the network data
  networks_start <- data.table::data.table(clusterL = network_modules$mergedColors, probeid = probe_id)
  
  # merge with affy info (contains entrez gene ids)
  affy_networks_merged <- merge(networks_start, affy, "probeid", all.x = TRUE) 
  
  # some probe ids have multiple gene symbols/entrez codes: gene1 /// gene2 /// ...
  
  # function to split the gene symbols/entrez code into different cells where each gene symbol gets its own row
  remove_dashes <- function(d){
    # determine the max number of genes per probe
    rounds <- max( str_count(d$entrez_gene, "///"), na.rm = TRUE )
    
    # generate a copy of the original data, change the column names for easy access in code
    d2 <- copy(d)
    setnames(d2, c("entrez_gene"), c("e1"))
    
    # splits the columns with multiple genes into columns with only one gene as many times as needed
    # output is a data frame that is wide with different genes 
    for(i in 1:rounds){
      d2 <- separate_(d2, paste0("e", i), paste0("e", i:(i+1)), sep = "///", extra = "merge")
    }
    
    # melt to get the original data with 1 gene_symbol and entrez column; may contain duplicate probe ids long
    d2 <- melt(d2, id.vars = c("probeid", "clusterL", "gene_symbol"), measure.vars = patterns("^e"), value.name = c("entrez_id"), na.rm = FALSE)
    
    # fix for clusters doing something funky with value name
    d2 <- dplyr::select(d2, probeid, clusterL, gene_symbol, matches("entrez_id"))
    setnames(d2, c("probeid", "clusterL", "gene_symbol", "entrez_id"))
    
    # formatting changes
    d2 <- d2 %>% 
      mutate(entrez_id = str_trim(entrez_id)) %>% 
      distinct() %>% 
      arrange(clusterL, gene_symbol)
    
    # return results
    return(d2)
  }
  
  # split probeid genes
  if( all(is.na(affy_networks_merged$entrez_gene)) ){
    
    out_data <- affy_networks_merged %>% 
      dplyr::select(probeid, clusterL, gene_symbol) %>% 
      distinct() %>% 
      arrange(clusterL, gene_symbol)
    
  } else{
    
    out_data <- remove_dashes(affy_networks_merged)
    
    # fix for NA entrez genes
    fix <- function(x){
      if( all(is.na(x$entrez_id)) ){
        return( x[1,] )
      } else{
        return( subset(x, !is.na(entrez_id)) )
      }
    }
    out_data <- out_data %>% 
      group_by(clusterL, gene_symbol) %>% 
      do(fix(.))    
  }
  
  # write to csv
  write.csv(out_data, file.path(output_folder, paste0(output_prefix,  "network.csv")), row.names = FALSE)
  
  # return results
  return(out_data)
  
}


# Run GO enrichment (preliminary)
#------------------------------------------------------------------------

# Function edit module eigengenes and saves results as csv
# Parameters: 
  # modules: output from cluster_modules
  # output_folder: folder for outputs
  # output_prefix: output file prefix
# Returns: a data frame

output_MEs <- function(modules, strains, output_folder = "", output_prefix = ""){
  
  # add strain names to MEs
  MEs <- modules$mergedMEs %>% 
    mutate(strain = strains) %>% 
    dplyr::select(strain, everything())
  
  # write to txt
  write.table(MEs, file.path(output_folder, paste0(output_prefix, "module_eigengenes.txt")), row.names = FALSE, sep = "\t", quote = FALSE)
  
  # return results
  return(MEs)
}

# Run GO enrichment (preliminary)
#------------------------------------------------------------------------

# Function to run GO enrichment analysis and saves results as csv
# Parameters: 
  # networkdf: output from network_to_df
  # output_folder: folder for outputs
  # output_prefix: output file prefix
# Returns: a data frame

runGO <- function(networkdf, output_folder = "", output_prefix = ""){
  
  # run GO enrichment analysis
  GOenr <- GOenrichmentAnalysis(networkdf$clusterL, networkdf$entrez_id, organism = 'mouse')
  
  # obtain best 10 terms in each module 
  enrichment <- GOenr$bestPTerms[[4]]$enrichment 
  
  # output results
  write.csv(enrichment, file.path(output_folder, paste0(output_prefix, "network_GOenrichment.csv")), row.names = FALSE)
  
  # return results
  return(enrichment)
}

