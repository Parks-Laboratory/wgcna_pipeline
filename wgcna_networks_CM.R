#------------------------------------------------------------------------
# Title: Run WGCNA for CM Data
# Date: July 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------


# CONTENTS:
  # Pull data from database and filter out genes
  # Run QC for missing data & outliers
  # Choose the soft-thresholding power - set params in separate wcgna_network_params.R
  # Run network (functions)
  # Format cluster data into data.frame
  # Save data

# NOTES:
  # Files to be transferred to clusters
    # wgcna pipeline files
    # affy_data.csv
    # network_expr_data.Rdata

#------------------------------------------------------------------------


# Input Command Line Args
#------------------------------------------------------------------------

library(optparse)

# set arguments
option_list <- list(
  make_option(c("--sql"), default = FALSE, help = "Whether the code is extracting data from SQL [default %default]"),
  make_option(c("--condor"), default = FALSE, help = "Whether code is running on Condor [default %default]"), 
  make_option(c("--i"), help = "What iteration to run, ranges from 0-3"), 
  make_option(c("--output_folder"), help = "Where to save outputs"),
  make_option(c("--opt_softpower"), default = FALSE, help = "Whether to run softpower analysis [default %default]")
)

# process arguments
opt <- parse_args(OptionParser(option_list = option_list))

# set input parameters
sql <- opt$sql
condor <- opt$condor
i <- as.numeric(opt$i) + 1
code_output_folder <- opt$output_folder
run_softpower <- opt$opt_softpower


# Load Libraries
#------------------------------------------------------------------------
library(reshape2)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(WGCNA)


# Run Networks
#------------------------------------------------------------------------

# run on condor for more memory options
if(sql){
  
  library(sqldf)
  library(RODBC)
  
  # load network code
  source("E:/John Li data/WGCNA/wgcna_network_pipeline.R")
  
  # set working directory
  setwd("E:/John Li data/WGCNA")
  
  
  # Pull Data from Database
  #------------------------------------------------------------------------
  
  # connect to SQL Server
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=external_labs;Trusted_Connection=Yes;DRIVER={SQL Server}')
  
  # function to extract CM data from database
  # Returns: list of data from sql for different tissues
  grab_sql <- function(){
    
    # obtain the names of table/views from SQL Server: tables that contain differences to run networks
    data_sources <- sqlTables(db) %>% 
      dplyr::select(TABLE_CAT, TABLE_SCHEM, TABLE_NAME) %>% 
      subset(str_detect(TABLE_NAME, "Chick_Munger_rna_expr"))
    
    # generate the SQL queries from the table names
    queries <- paste("select mouse_id as strain, probe_id, norm_expression as expression from", paste(data_sources$TABLE_SCHEM, data_sources$TABLE_NAME, sep = "."), "order by probe_id, strain")
    
    # connect to SQL Server and extract data
    data <- lapply(queries, function(q) sqlQuery(db, q))
    names(data) <- str_extract(queries, "(chow|HF)_.")
    
    # return sql results
    return(data)
  }
  
  # extract data
  sql_data <- grab_sql()
  
  # close the database
  odbcClose(db)
  
  
  # Process SQL Data
  #------------------------------------------------------------------------
  
  # process SQL data
  process_sql_data(sql_data = sql_data, use_col = "expression", outlier_cutoff = 200, output_folder = code_output_folder)
  
  
  # Choose the soft-thresholding power
  #------------------------------------------------------------------------
  
  # load data
  load( file.path(code_output_folder, "network_expr_data.Rdata") ) # called expr_data and strain_names
  options <- names(expr_data)
  
  # pick_softpower
  if(run_softpower){
    
    for(opt in options){
      pick_soft_power(expr_data[[opt]], output_folder = code_output_folder, output_suffix = paste0("_", opt))
    }
    
  }
  
  
  # Obtain Affy Data
  #------------------------------------------------------------------------
  
  # import the affy annotation data
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=external_labs;Trusted_Connection=Yes;DRIVER={SQL Server}')
  affy <- sqlQuery(db, "select ensembl_id as probeid, gene_symbol, entrez_gene from dbo.Chick_Munger_rna_probe_annotations")
  odbcClose(db)
  
  write.csv(affy, file.path(code_output_folder, "affy_data.csv"), row.names = FALSE)
  
  
} else if(condor){
  
  
  # Load Data
  #------------------------------------------------------------------------
  
  # load network code
  source("wgcna_functions_condor.R")
  source("wgcna_network_pipeline.R")
  
  # load gene annotations
  affy <- fread("affy_data.csv")
  
  # load data
  load("network_expr_data.Rdata") # expr_data & strain_names
  options <- names(expr_data)
  
  # set params
  kind <- options[i]
  expr_data <- expr_data[[kind]]
  strain_names <- strain_names[[kind]]
  
  
  # Run Pipeline
  #------------------------------------------------------------------------

  # set soft power values
  sp <- switch(
    kind,
    "chow_M" = 4,
    "chow_F" = 2,
    "HF_M" = 3,
    "HF_F" = 3
  )
  
  run_wgcna(
    data_id = kind, 
    expr_data = expr_data, 
    probeid = colnames(expr_data), 
    strains = strain_names, 
    affy = affy, 
    # clinical_traits = , NO CLINICAL TRAITS
    
    soft_power = sp,
    min_mod_size = 30,
    deep_split = 3,
    ME_cor_thres = 0.1,
  
    output_folder = code_output_folder,
    output_genetree = TRUE
  )

  
# second portion of code just makes plots
} else{
  
  # load network code
  source("E:/John Li data/WGCNA/wgcna_network_pipeline.R")
  
  # set working directory
  setwd("E:/John Li data/WGCNA")
  
  # load the data: expression data and Rdata from Condor
  load(file.path(code_output_folder, "network_expr_data.Rdata"))
  
  options <- names(expr_data)
  
  for(kind in options){
    
    # load data
    load(file.path(code_output_folder, paste0(kind, "_networkConstruction.Rdata")))
    load(file.path(code_output_folder, paste0(kind, "_gene_module_data.Rdata")))
    
    # print merged modules
    png(file.path(code_output_folder, paste0("clustering_after_merge_", kind, ".png")), height = 500, width = 700)
    plotDendroAndColors(genes$geneTree, cbind(genes$dynamicColors, modules$mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    dev.off()
     
  }
  
}

