#------------------------------------------------------------------------
# Title: Template for Running WGCNA on Clusters
# Date: May 17 2016
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

#------------------------------------------------------------------------


# Input Command Line Args - TODO Add as Necessary
#------------------------------------------------------------------------

library(optparse)

# set arguments
option_list <- list(
  make_option(c("--sql"), default = FALSE, help = "Whether to run SQL to obtain data (only need to run once) [default %default]"),
  make_option(c("--output_folder"), help = "Folder for code outputs"), 
  make_option(c("--opt_softpower"), default = FALSE, help = "Whether to run softpower analysis [default %default]")
)

# process arguments
opt <- parse_args(OptionParser(option_list = option_list))

# set input parameters
sql <- opt$sql
run_softpower <- opt$opt_softpower


# TODO: set the working directory here
  # the working directory contains the R scripts
  # don't forget to create a (separate) folder for code outputs
setwd("E:/MICROARRAY DATA/CHANGE/") 
code_output_folder <- opt$output_folder
if(!dir.exists(code_output_folder)) dir.create(code_output_folder)


# Load Libraries
#------------------------------------------------------------------------
library(reshape2)
library(stringr)
library(plyr)
library(dplyr)
library(dtplyr)
library(tidyr)
library(sqldf)
library(data.table)
library(RODBC)
library(WGCNA)
source("E:/MICROARRAY DATA/wgcna_pipeline/wgcna_network_pipeline.R")


# Process SQL Data - TODO: edit this process manually to select and filter the data from SQL
#------------------------------------------------------------------------

if(sql){
  
  # Pull Data from Database
  #------------------------------------------------------------------------
  
  # connect to SQL Server
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
  
  # Returns: list of data from sql for different tissues
  grab_sql <- function(){
    
    # obtain the names of table/views from SQL Server: tables that contain differences to run networks
    data_sources <- sqlTables(db) %>% 
      subset(TABLE_SCHEM == "dbo") %>% 
      dplyr::select(TABLE_CAT, TABLE_SCHEM, TABLE_NAME) %>% 
      subset(str_detect(TABLE_NAME, "str_match")) # TODO: edit str_match to match names of table in SQL server
    
    # generate the SQL queries from the table names TODO: enter the column names here
    queries <- paste("select strain_id as strain, probe_id as probe_id, expression_values as expression from", paste(data_sources$TABLE_SCHEM, data_sources$TABLE_NAME, sep = "."), "order by probe_id, strain")
    
    # connect to SQL Server and extract data
    data <- lapply(queries, function(q) sqlQuery(db, q))
    names(data) <- str_extract(data_sources$TABLE_NAME, "") # TODO: extract names of data entries; usually is related to the table_name
    
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
  
}


# Load Data
#------------------------------------------------------------------------

load( file.path(code_output_folder, "network_expr_data.Rdata") ) # called expr_data and strain_names


# Choose the soft-thresholding power
#------------------------------------------------------------------------

options <- names(expr_data)
if(run_softpower){
  
  for(opt in options){
    pick_soft_power(expr_data[[opt]], output_folder = code_output_folder, output_suffix = paste0("_", opt))
  }

}


# Run network (code) TODO: edit data as necessary
#------------------------------------------------------------------------

if(!run_softpower){
  
  # import the affy annotation data TODO: edit affy data as needed 
    # this should contain columns probeid, gene_symbol, entrez_gene for GO enrichments
    # ok if there is no entrez_gene column available; set entrez_gene to be NA or (in sql statement) NULL as entrez_gene
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
  affy <- sqlQuery(db, "select probesetID as probeid, gene_symbol, entrez_gene from dbo.affy_microarray_annotation") 
  odbcClose(db)
  
  # import clinical trait differences TODO: edit clinical trait data if needed
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
  clinical_traits <- sqlQuery(db, "select * from clinical_trait_differences_male") 
  odbcClose(db)
  
  # run pipeline on all kinds
  for(opt in options){
    
    # TODO: set soft_power based on outputs for all options
    sp <- 4
    
    run_wgcna(
      data_id = opt, 
      expr_data = expr_data[[opt]], 
      probeid = colnames(expr_data[[opt]]), 
      strains = strain_names[[opt]], 
      affy = affy,  
      clinical_traits = clinical_traits, # TODO: if don't want clinical trait/module eigengene cor or no clinical trait data, comment out this line
      soft_power = sp, 
      min_mod_size = 30, 
      deep_split = 3, 
      ME_cor_thres = 0.1, 
      output_folder = code_output_folder,
      output_genetree = FALSE
    )
    
  }
}

