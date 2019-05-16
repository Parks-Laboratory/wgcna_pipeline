# Generate R tarball for upload to condor
* Go to Interactive Condor and do the following:

1. [John@build]$ tar -xzf R-3.6.0.tar.gz
2. [John@build]$ cd R-3.6.0
3. [John@build]$ ./configure --prefix=$(pwd)
4. [John@build]$ make
5. [John@build]$ make install
6. [John@build]$ cd ..

* Open R program in condor:

7. [John@build]$ R-3.6.0/lib64/R/bin/R 
 > install.packages('package_name')   ...... 
   (needed packages names in library.txt file, RODBC is not availabe since it needs sudo authoriaty)
 
 exit the R
 > q()   
  
* change RHOME directory
8. [John@build]$ nano R-3.6.0/lib64/R/bin/R

change from : R_HOME_DIR=/var/lib/condor/execute/slot1/dir_554715/R-3.1.0/lib64/R
change to: R_HOME_DIR=$(pwd)/R

9.  [John@build]$ mv R-3.6.0/lib64/R ./
10. [John@build]$ tar -czvf R.tar.gz R/
11. [John@build]$ exit 

* The R.tar.gz will be transfered back to home directory after exit condor
* R.tar.gz (with all necessary libraries for WGCNA) is in E:/condor_tools


# Running the WGCNA Pipeline

For example of running on our server: see REACTIVE_HDMP_NETWORKS  
For example of running on clusters: see INDIVIDUAL_NETWORKS or ChickMungeretal2016

If the number of probes is large > 4000, then need to run code separately on the clusters

When running pipeline on the clusters, transfer the following files
* corresponding network_expr_data.Rdata
* corresponding affy_data.csv
* all files in wgcna pipeline folder
* wgcna.sh
* wgcna.sub

On the cluster  

1. make a separate directory called **wgcna/**
2. make a separate directory called **node_files/** inside the **wgcna/** directory  
3. transfer all R, Rdata, csv files into **node_files/** directory  
4. edit the **wgcna.sh** file: in the last line change the name of the wgcna_network script to the correct file  
5. edit the **wgcna.sub** file: change the queue number to the number of data sets (ie the length of the list expr_data, contained inside network_expr_data)  
6. copy **R.tar.gz** from cluster into the **wgcna/** directory  
7. run the command `condor_submit wgcna.sub`
