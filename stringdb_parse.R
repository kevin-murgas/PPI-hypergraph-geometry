# script for parsing STRINGdb PPI data
# load in PPI using STRINGdb package
# format and save as node + edge lists (for MATLAB topology script)

library(STRINGdb)
library(tidyverse)

# set working directory as needed, where to save the string node/edge files
# setwd(your_working_directory)
library(here)

# load STRING database
string_db <- STRINGdb$new(version="11", species=9606,
                          score_threshold=900, input_directory="")

# get protein names and PPI graph
proteins <- string_db$get_proteins()
g <- string_db$get_graph() %>% simplify

library(igraph)
gorder(g) # 12396 nodes
gsize(g) # 324152 edges (undirected)

# take largest strongly connected component
components = clusters(g,"strong")
max_comp_id = which.max(components$csize)
vind = V(g)[components$membership == max_comp_id]
g_max = induced_subgraph(g,vind)
proteins_max = proteins %>% filter(protein_external_id %in% V(g_max)$name)

gorder(g_max) # 11919 nodes
gsize(g_max) # 323823 edges

# convert string_id to gene names
gene_names = proteins_max$preferred_name

# save nodes and edges as csv
gene_names_df <- data.frame(gene=gene_names) # could add annotations, alternate gene names
write_csv(gene_names_df,"stringdb11_nodes.csv")

# get edge list and reorder
edge_df = as_edgelist(g_max, names=F) %>% as.data.frame() %>%
  mutate(V1 = gene_names[V1], V2 = gene_names[V2]) %>%
  rename(V3 = V2) %>% add_column(V2="interacts_with", .before=2)
write_csv(edge_df,"stringdb11_edges.csv")


## SPARSIFICATION
  # Banerji 2013 sparsification procedure:
  # they used this procedure to annotate genes as EC MR or IC, favor "external domain"
  # extracellular space/region using the following set of GO-terms:
  # GO:0005102, GO:0008083,GO:0005125, GO:0005615, GO:0005576, GO:0044421
  #  transmembrane receptor activity (GO:0004888, GO:0004872, GO:0005886)
  # intracellular domain or to biological functions associated with the intracellular region
  # GO:0005622, GO:0044424, GO:0005634, GO:0005737, GO:0005829, GO:0000139, GO:0035556, GO:0007243, GO:0006468)
  # only allow EC-EC, EC-MR/MR-EC, MR-IC/IC-MR, IC-IC
  # then take maximally connected component

# annotate proteins with GO location, using mygene (query may take a minute)
library(mygene)
res <- queryMany(gene_names, scopes='symbol', fields=c('go'), species='human') %>% as.data.frame

# summarize locations to extracellular, membrane, or intracellular
# only use annotations from experimental evidence
# EC, MR, IC tags defined in Banerji 2013
evidence_tags = c("EXP","IDA","IPI","IMP","IGI","IEP","HTP","HDA","HMP","HGI","HEP")
ec_tags = c("GO:0005102","GO:0008083","GO:0005125", "GO:0005615", "GO:0005576", "GO:0044421")
mr_tags = c("GO:0004888","GO:0004872","GO:0005886")
ic_tags = c("GO:0005622","GO:0044424","GO:0005634", "GO:0005737", "GO:0005829", "GO:0000139", "GO:0035556", "GO:0007243", "GO:0006468")

# summarize subcellular location for each gene
cc_term = gene_names
for (i in 1:length(gene_names)) { # loop may take a few minutes
  if (!(i %% 1000)) { print(sprintf("%i/%i",i, length(gene_names))) }
  gocc = res %>% dplyr::filter(query==gene_names[i]) %>% pull(go.CC)
  l = sapply(gocc,length) # detect length in case of duplicate entries (usually only 1 has go annotations)
  if (!sum(l)) { cc_term[i] = "na"; next } # use na if not annotated
  gocc = gocc[[which(l>0)]] %>% as.data.frame
  gocc = gocc %>% dplyr::filter(evidence %in% evidence_tags) # only use experimental evidence annotations
  term = "ua" # initially unknown ua
  if (any(gocc$id %in% ec_tags)) {term = "ec"}
  else if (any(gocc$id %in% mr_tags)) {term = "mr"}
  else if (any(gocc$id %in% ic_tags)) {term = "ic"}
  cc_term[i] = term
}
cc_term %>% as.factor %>% levels   # "ec" "ic" "mr" "na" "ua"
cc_term %>% as.factor %>% tabulate #  631 3082 1288  708 6179

# remove edges if they are ec-ic or ic-ec, otherwise leave alone
g_sparse <- delete.edges(g_max, E(g_max)[ V(g_max)[which(cc_term == 'ec')] %--%
                                            V(g_max)[which(cc_term == 'ic')] ])
g_sparse <- delete.edges(g_sparse, E(g_sparse)[ V(g_sparse)[which(cc_term == 'ic')] %--%
                                                  V(g_sparse)[which(cc_term == 'ec')] ])
adj = as_adjacency_matrix(g_sparse)
adj[which(cc_term=="ec"), which(cc_term=="ic")] %>% sum() # check is zero
adj[which(cc_term=="ic"), which(cc_term=="ec")] %>% sum()

# sparsification complete
# take max connected component
components = clusters(g_sparse,"strong")
max_comp_id = which.max(components$csize)
vind = V(g_sparse)[components$membership == max_comp_id]
g_sparse_max = induced_subgraph(g_sparse,vind)
proteins_max = proteins %>% filter(protein_external_id %in% V(g_sparse_max)$name)
gene_names = proteins_max$preferred_name

gorder(g_sparse_max) # 11888 nodes
gsize(g_sparse_max) # 315130 edges

# save nodes and edges as csv
gene_names_df <- data.frame(gene=gene_names)
write_csv(gene_names_df,"stringdb11_sparse_nodes.csv")

edge_df = as_edgelist(g_sparse_max, names=F) %>% as.data.frame() %>%
  mutate(V1 = gene_names[V1], V2 = gene_names[V2]) %>%
  dplyr::rename(V3 = V2) %>% add_column(V2="interacts_with", .before=2)
write_csv(edge_df,"stringdb11_sparse_edges.csv")


# to load sparsified STRINGdb
#gene_names = read_csv("stringdb11_sparse_nodes.csv") %>% pull(gene)
#edge_df = read_csv("stringdb11_sparse_edges.csv")
