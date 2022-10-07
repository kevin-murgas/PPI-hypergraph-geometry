# PPI-hypergraph-geometry
Code for publication:\
Hypergraph Geometry Reflects Higher-Order Dynamics in Protein Interaction Networks\
Kevin A. Murgas, Emil Saucan, Romeil Sandhu\
Scientific Reports (under review)

## Requirements:
- MATLAB
- R

## Summary of code script files:
1. **stringdb_parse.R**, **convert2mat.m** - These two scripts prepare the PPI topology from the STRINGdb protein-protein interaction database. The R script accesses the database and saves node- and edge-list files; the MATLAB script then processes the node- and edge-list into an adjacency matrix and list of gene names (as HGNC gene symbol), which are saved as a single .mat file containing a standardized representation of a PPI topology.
2. **examine_topology.m** - This MATLAB script is used to analyze topological aspects of the 1-dimensional (1D) graph model and 2-dimensional (2D) hypergraph model of the PPI network, including numbers of vertices, edges, and faces (in the 2D case), as well as topological characteristics including the Euler number. Subsequently, edge- and face-degree distributions are assessed.
3. **preprocess_expression.m** - This MATLAB script preprocesses gene expression data (from various RNA-seq datasets from public databases such as GEO, TCGA, etc) and formats the normalized expression into a standardized .mat file containing the expression matrix along with sample (column) and gene (row) names. The normalization procedure includes quantile normalization and log_2 transform. Because of the variability in datasets, this processing script incorporates user input to ensure the correct data columns are selected and optional proecessing steps can be omitted. Note, because of the challenge of different data formats, some user tweaking of the code may be required to correctly prepare their own datasets.
4. **curvature_script.m** - This MATLAB function script loads in specified topology and expression .mat files (prepared using the above scripts 1 and 3) corresponding to a PPI topology to build the network model and a gene expression dataset to be computed for geometric network analysis. The algorithm first combines the data by intersecting genes shared in both the gene expression data and PPI topology, then taking the maximally connected network component. Then, the algorithm runs through each sample (column) of gene expression, constructing a weighted network model (see manuscript for details) and then computing geometric features of the model, including Forman-Ricci curvature of the 1D graph model and 2D hypergraph extended model, as well as graph entropy. Scalar contraction on the vertices and a global average are computed. Finally these results from all samples are saved in .csv format.
   * This code is optimized with parallel processing to compute FR-curvature the edges of the network in the 2D model. For this reason we highly suggest the use of a high-performance compute cluster. Because the MATLAB is a function script, input arguments are expected when calling from BASH commands; below is an example of how to call the script:
      ```console
      matlab -batch "curvature_script $BASEDIR $EXPRFILE $TOPOFILE $OUTDIR"
      ```
   * where $BASEDIR, $EXPRFILE, $TOPOFILE, $OUTDIR specify respectively the base directory (where to run code from), expression file (.mat, prepared by preprocess_expression.m), topology file (prepared by convert2mat.m), and output directory to save results files.
5. **analyze2Dresults.R** - This R script is used to examine the results from the curvature_script.m analysis. The script loads in results files from a specified experiment (i.e. gene expression dataset). Then various data are represented in figures corresponding to the figures in our manuscript.

## Support
Any questions about implementation or code bugs are welcome. Please use Github issues, or email the authors.
