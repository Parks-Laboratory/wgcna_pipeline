# Running the WGCNA Pipeline

For example of running on our server: see REACTIVE_HDMP_NETWORKS

For example of running on clusters: see INDIVIDUAL_NETWORKS or ChickMungeretal2016

If the number of probes is large > 4000, then need to run code separately on the clusters

When running pipeline on the clusters, transfer the following files
- corresponding network_expr_data.Rdata
- corresponding affy_data.csv
- all files in wgcna pipeline folder
- wgcna.sh
- wgcna.sub

On the cluster
1) make a separate directory called wgcna/
2) make a separate directory called node_files/ inside the wgcna/ directory
3) transfer all R, Rdata, csv files into node_files/ directory
4) edit the wgcna.sh file: in the last line change the name of the wgcna_network script to the correct file
5) edit the wgcna.sub file: change the queue number to the number of data sets (ie the length of the list expr_data, contained inside network_expr_data)
6) copy R.tar.gz from cluster into the wgcna/ directory
7) run the command **condor_submit wgcna.sub**
