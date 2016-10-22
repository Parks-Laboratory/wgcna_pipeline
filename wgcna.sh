#!/bin/bash

# untar your R installation
tar -xzf R.tar.gz

# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH

mkdir outputs

# run R, with the name of your  R script
R CMD BATCH '--args --sql=FALSE --condor=TRUE --i='$1' --output_folder=outputs' wgcna_networks_CM.R wgcna_output_$1.Rout

mv outputs/* .
