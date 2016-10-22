#------------------------------------------------------------------------
# Title: Run WGCNA for Individual Expression Data
# Date: May 23 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------


# CONTENTS:
  # Input Command Line Args
  # Process SQL Data
  # Pick Soft Power
  # Run WGCNA

# NOTES:
  # Search TODO
  # Files to be transferred to clusters
  # wgcna pipeline files
  # affy_data.csv 
  # clinical_traits.csv (if necessary)
  # network_expr_data.Rdata

#------------------------------------------------------------------------


# Input Command Line Args TODO: adjust as necessary
#------------------------------------------------------------------------

library(optparse)

# set arguments
option_list <- list(
  make_option(c("--sql"), default = FALSE, help = "Whether the code is extracting data from SQL [default %default]"),
  make_option(c("--condor"), default = FALSE, help = "Whether code is running on Condor [default %default]"), 
  make_option(c("--i"), help = "What iteration to run, ranges from 0-5"), 
  make_option(c("--output_folder"), help = "Where to save output"),
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
  source("E:/MICROARRAY DATA/wgcna_pipeline/wgcna_network_pipeline.R")
  
  # set working directory TODO: set working directory
  setwd("E:/MICROARRAY DATA/CHANGE/")
  
  
  # Pull Data from Database
  #------------------------------------------------------------------------
  
  # connect to SQL Server
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
  
  # function to extract highfat v chow expression difference data
  # Returns: list of data from sql for different tissues
  grab_sql <- function(){
    
    # obtain the names of table/views from SQL Server: tables that contain differences to run networks
    data_sources <- sqlTables(db) %>% 
      subset(TABLE_SCHEM == "dbo") %>% 
      dplyr::select(TABLE_CAT, TABLE_SCHEM, TABLE_NAME) %>% 
      subset(str_detect(TABLE_NAME, "avg_expression_by_strain")) %>% # TODO: edit str_match to match names of table in SQL server

    # generate the SQL queries from the table names TODO: enter the column names here
    queries <- paste("select strain_id as strain, probe_id as probe_id, expression_values as expression from", paste(data_sources$TABLE_SCHEM, data_sources$TABLE_NAME, sep = "."), "order by probe_id, strain")
      
    # connect to SQL Server and extract data
    data <- lapply(queries, function(q) sqlQuery(db, q))
    names(data) <- data_sources$TABLE_NAME # TODO: extract names of data entries; usually is related to the table_name
    
    # return sql results
    return(data)
  }
  
  # extract data
  sql_data <- grab_sql()
  
  # close the database
  odbcClose(db)
  
  # TODO: if more processing of data is required, do it here
  
  
  # Process SQL Data
  #------------------------------------------------------------------------
  
  # process SQL data 
    # TODO: edit the use_col based on the column name representing expression in the data set (see the query)
    # TODO: edit outlier_cutoff based on outlier hierchical clustering plots; run through first with 100 - it is usually ok, but may change if needed
  process_sql_data(sql_data = sql_data, use_col = "expression", outlier_cutoff = 100, output_folder = code_output_folder)

  
  # Softpower Analysis
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
  
  # import the affy annotation data TODO: edit affy data as needed 
    # this should contain columns probeid, gene_symbol, entrez_gene for GO enrichments
    # ok if there is no entrez_gene column available; set entrez_gene to be NA or (in sql statement) NULL as entrez_gene
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
  affy <- sqlQuery(db, "select probesetID as probeid, gene_symbol, entrez_gene from dbo.affy_microarray_annotation")
  odbcClose(db)
  
  write.csv(affy, file.path(code_output_folder, "affy_data.csv"), row.names = FALSE)
  
  
  # Obtain Clinical Traits
  #------------------------------------------------------------------------
  
  # import clinical trait differences TODO: edit clinical trait data if needed
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
  clinical_traits <- sqlQuery(db, "select * from clinical_trait_differences_male") 
  odbcClose(db)
  
  write.csv(clinical_traits, file.path(code_output_folder, "clinical_trait_data.csv"), row.names = FALSE)
  
  
} else if(condor){
  
  # Load Data
  #------------------------------------------------------------------------
  
  # load network code
  source("wgcna_functions_condor.R")
  source("wgcna_network_pipeline.R")
  
  # load gene annotations TODO: if data not available, comment out these lines
  affy <- fread("affy_data.csv")
  clinical_traits <- fread("clinical_trait_data.csv")
  
  # load data
  kind <- options[i]
  load("network_expr_data.Rdata") # expr_data & strain_names
  options <- names(expr_data)
  
  # set params
  kind <- options[i]
  expr_data <- expr_data[[kind]]
  strain_names <- strain_names[[kind]]
  
  
  # Run Pipeline
  #------------------------------------------------------------------------
  
  # TODO: set soft_power based on outputs
  sp <- 4
  
  run_wgcna(
    data_id = kind, 
    expr_data = expr_data, 
    probeid = colnames(expr_data), 
    strains = strain_names, 
    affy = affy,  
    clinical_traits = clinical_traits, # TODO: if don't want clinical trait/module eigengene cor or no clinical trait data, comment out this line
    
    soft_power = sp,
    min_mod_size = 30,
    deep_split = 3,
    ME_cor_thres = 0.1,
    
    output_folder = code_output_folder,
    output_genetree = TRUE
  )
  
  
  # second portion of code requires db access - run on server
} else{
  
  # load network code
  source("E:/MICROARRAY DATA/wgcna_pipeline/wgcna_network_pipeline.R")
  
  # set working directory TODO set working directory
  setwd("E:/MICROARRAY DATA/CHANGE/")
  
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

