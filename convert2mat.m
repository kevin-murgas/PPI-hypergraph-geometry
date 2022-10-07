% convert PPI network topology data into MATLAB object
% adjacency matrix and corresponding gene names
% save output containing:
% adj: nxn adjacency matrix (double, but only 0/1)
% names: nx1 list of gene names

% load in topology
% Using STRINGdb v11 (sparsified)
% node and edge list prepared with stringdb_parse.R script
topodir = "STRINGdb_data";
nodefile = "stringdb11_sparse_nodes.csv";
edgefile = "stringdb11_sparse_edges.csv";

fprintf('Load node and edge CSV files\n');
nodes = readtable(fullfile(topodir,nodefile));
nodes = string(table2array(nodes)); % convert to string array
edges = readtable(fullfile(topodir,edgefile));

fprintf('Match to node indices, build adjacency matrix\n')
[ml,left_nodes] = ismember(edges.V1,nodes); % match left side of edge gene names (column V1) to gene names in nodes
[mr,right_nodes] = ismember(edges.V3,nodes); % again for right side of edges (V3)
mb = and(ml,mr); % edges with both left and right nodes present in expression data genes
left_nodes = left_nodes(mb); right_nodes = right_nodes(mb);
adj_all = sparse(left_nodes, right_nodes, 1, length(nodes), length(nodes));
adj_all = adj_all + adj_all'; % make symmetric (for STRINGdb, undirected)

fprintf('Take maximally connected component\n')
G = graph(adj_all);
[connind, binsize] = conncomp(G); % labels all nodes by connected component, and gives # nodes per component
maxconnind = find(connind==find(binsize==max(binsize))); % find inds where connectivity label is that of the largest bin size
adj = adj_all(maxconnind,maxconnind);
topoNames = nodes(maxconnind);

% save adj and names in one .mat file
fprintf('Saving .mat file\n')
savedir = "Topology";
filename = "topo_stringdb11sparse.mat";
save(fullfile(savedir,filename),'adj','topoNames')