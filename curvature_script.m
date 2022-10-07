% Forman Curvature Analysis Script
% By Kevin Murgas, Stony Brook Univ. Biomedical Informatics

% takes in simple command line arguments
% input expression data: single file (preprocessed)
% topology: .mat file of PPI vertex names and adjacency matrix
% output: folder directory

% compute entropy, forman curvature in 1D and 2D

function curvature_script(base_dir, expr_file, topo_file, output_dir)
%% read args in from command line

% manually defined defaults
if nargin<1
    base_dir = pwd;
end
if nargin<2
    expr_file = 'Data_proc/GSE75748_3_ql2.mat';
end
if nargin<3
    topo_file = 'Topology/topo_stringdb11sparse.mat';
end
if nargin<4
    output_dir = 'Analysis/job';
end
fprintf("\nNetwork Curvature Script - Kevin Murgas, Stony Brook University Dept. of Biomedical Informatics\nbase_dir = %s\nexpr_file = %s\ntopo_file = %s\noutput_dir = %s\n",base_dir,expr_file,topo_file,output_dir)

cd(base_dir);

%% load data
% load in the expression and topology data files
% inner join to produce variables: adj_max, genes_max, expr_max

% load expression data (pre-processed RNA expression data, contains exprNorm, geneName, sampName)
load(expr_file);
nSamp = length(sampName);

% load topology (as .mat containing adj and topoNames)
load(topo_file)

% inner join expression and topology genes
[~, idx] = ismember(intersect(geneName, topoNames), geneName); % get indices of gene names in both expression and topology data
exprNames = geneName(idx);
expr_join = exprNorm(idx,:);
[~, idx] = ismember(exprNames, topoNames);% take only rows of adjacency matrix corresponding to vertices
adj_all = adj(idx,idx);

% take maximal connected component indices from network
G = digraph(adj_all);
[connind, binsize] = conncomp(G); % labels all vertices by connected component, and gives # vertices per component
maxconnind = find(connind==find(binsize==max(binsize))); % find inds where connectivity label is that of the largest bin size

% adj_max, expr_max are the combined PPI network and gene expression
adj_max = adj_all(maxconnind,maxconnind);
genes_max = exprNames(maxconnind);
expr_max = expr_join(maxconnind,:);

fprintf("Matched expression data to PPI. Size of maximum connected component = %i x %i\n",size(adj_max))

%% prepare for processing
% compute some topology-dependent variables ahead of curvature computation

adj = logical(adj_max); % use maximally-connected adjacency matrix
nv = size(adj,1);
ne = sum(adj(:));

%deg_in = sum(adj,1);
deg_out = sum(adj,2);

% compute left/right vertex pairs from adj
[lv, rv] = ind2sub(size(adj), find(adj));

% index edges and reverse edge pairs (topology-dependent)
% ei: edge index in adjacency matrix, er: reverse edge index
% rev: binary indicator if opposite edge is present
% revi: pointer to which edge is the reverse pair
fprintf('Index edges and reverse edge pairs\n')
ei = sub2ind([nv nv], lv, rv); % ind of each edge in adj
er = sub2ind([nv nv], rv, lv); % ind of corresponding reverse edge
[rev, revi] = ismember(ei,er); % find which positions have a matching reverse index
revi = revi+~rev; % add 1 for the spaces with 0 to make work as index array

% directed graph faces
% face identification:
% make an edge-list of faces                                to visualize directed faces:  > v3 \                    / v3 <     :
% each edge (lv/rv pair) stores a third vertex                                           /      \                  /      \    :
% to represent (positive or negative) directed faces                                    /        V                V        \   :
% automatically encode + or - face in two different lists                    positive: v1 ------> v2  negative: v1 ------> v2  :
% ideal to summarize and vectorize for edge curvature                                     edge i                   edge i      :

% get neighboring vertices
nbrs_in = accumarray(rv, lv,[nv 1], @(x) {x}); % accumulate cell array of ingoing neighbors at each vertex
nbrs_out = accumarray(lv, rv,[nv 1], @(x) {x}); % accumulate cell array of outgoing neighbors at each vertex

% then use the intersect function to find face vertices at each edge
fprintf('Get intersecting triangular face vertices\n')
f_pos_cell = cell(1,length(lv));
f_neg_cell = f_pos_cell;
tic
for ie = 1:length(lv)
    % get end vertices of edge
    v1 = lv(ie); v2 = rv(ie); % get vertices of edge
    
    % find incident faces by identifying vertices that form triangles with edge ie
    f_pos_cell{ie} = intersect(nbrs_out{v1}, nbrs_in{v2}); % positive directed faces
    f_neg_cell{ie} = intersect(nbrs_in{v1}, nbrs_out{v2}); % negative directed
end
toc
fprintf('%i + faces, %i - faces identified\n', sum(cellfun(@length, f_pos_cell)), sum(cellfun(@length, f_neg_cell)))


%% process each sample
% with use of parallel processing for 2D curvature of edges

% create parallel pool
delete(gcp('nocreate')) % close any current pool
numCores = feature('numcores');
warning('off','MATLAB:datetime:NonstandardSystemTimeZone'); % this line is to suppress an annoying timezone warning
parpool(numCores); % new pool with all available cores

fprintf('Begin curvature analysis (%i samples)\n',nSamp)
sampEntropy = struct;
sampCurvature_1D = struct;
sampCurvature_2D = struct;

t1=tic;
for iSamp = 1:nSamp
    fprintf('\nSample %i/%i\n',iSamp,nSamp)
    ti = tic;
    
    % define vertex+edge weights
    v_wts = expr_max(:,iSamp);
    e_wts = (v_wts(lv).*v_wts(rv)); % mass action
    e_wt_mat = sparse(lv, rv, e_wts, nv,nv);
    sto_mat = e_wt_mat./sum(e_wt_mat,2); % normalize outgoing edge weights
    
    % local network entropy
    Slocal = sum(spfun(@(x) -x.*log(x), sto_mat),2);
    
    % define geometric edge weights
    e_wts = 1./(deg_out(lv).*sto_mat(ei)); % resistance metric (Gobel & Jagers 1974)
    e_wt_mat = sparse(lv, rv, e_wts, nv,nv);
    
    % precompute vertex/sqrt(edge) sums for each vertex, sample-dependent
    invsqrtedge_mat = sparse(lv, rv, 1./sqrt(e_wts), nv,nv); % reshape 1/sqrt(edge weights) into adjacency matrix
    vsum_in = v_wts.*sum(invsqrtedge_mat,1)'; % take sum over end vertices (inward)
    vsum_out = v_wts.*sum(invsqrtedge_mat,2); % take sum over start vertices (outward)
    
    % network curvature
    % Forman-Ricci 1D formula (directed)
    Clocal_1D = v_wts(lv) + v_wts(rv) - sqrt(e_wts).*(vsum_in(lv) + vsum_out(rv) - rev.*(v_wts(lv) + v_wts(rv))./sqrt(e_wts(revi)));
    
    % Forman-Ricci 2D formula
    Clocal_2D = 0*Clocal_1D;
    parfor ie = 1:ne
        we = e_wts(ie); % weight of given edge
        
        % get end vertices of edge
        v1 = lv(ie); v2 = rv(ie); % get vertices of edge
        wv = v_wts([v1 v2])';
        
        % find incident faces by identifying vertices that form triangles with edge ie
        f_pos = f_pos_cell{ie}; % positive directed faces (3rd vertex)
        f_neg = f_neg_cell{ie}; % negative directed
        
        % face weights by mass action
        % w(f) = ( w(e1) * w(e2) * w(e3) ) ^ (2/3) = (geometric mean of edge weights)^2
        % use indices of lv, rv, and 3rd vertex in f_pos_cell{ie}
        
        % geometric mean squared of edge weights
        wf = [geomean([repmat(e_wt_mat(v1,v2), size(f_pos)) e_wt_mat(v1,f_pos)' e_wt_mat(f_pos,v2)],2);...
              -geomean([repmat(e_wt_mat(v1,v2), size(f_neg)) e_wt_mat(v2,f_neg)' e_wt_mat(f_neg,v1)],2)].^2;

        % sum term 1 = edge/face ratio
        sum1 = sum(we./wf,'all');
        % sum term 2 = vertex/edge ratio
        sum2 = sum(wv./we);
        
        % (-) sum term 3: parallel components = edges sharing one vertex or one face but not both
        % set A shares a face but no vertex
        % by definition of faces as triangles sharing edge ie, there are no edges that share a face but not one of the two vertices
        sum3a = 0; % for triangle faces, no edges are parallel via face
        
        % set B shares a vertex but no face
        % edges that are concordant in direction (into n1 or out of n2) are positive, opposite direction edges are negative
        % compute using vsum_in/out then removing unneeded edges
        % exclude: ei, f_neg->n1, n2->f_neg, and reverse edge if present
        % subtract these terms out, as they are included in vsum
        sum3b = vsum_in(v1) + vsum_out(v2) - sum(v_wts(v1)./sqrt(e_wt_mat(f_neg,v1))) - sum(v_wts(v2)./sqrt(e_wt_mat(v2,f_neg)));
        if adj(v2,v1)
            sum3b = sum3b - sum(wv)/sqrt(e_wt_mat(v2,v1));
        end
        sum3b = sum3b / sqrt(we);
        
        % put it all together into Forman curvature
        Clocal_2D(ie) = we*(sum1 + sum2 - (sum3a + sum3b));
    end
    
    % convert to  vertex curvature (sum of all edges in/out of vertex)
    Cv_out_1D = accumarray(lv, Clocal_1D);
    Cv_in_1D = accumarray(rv, Clocal_1D);
    %Cv_1D = Cv_in_1D - Cv_out_1D;
    Cv_out_2D = accumarray(lv, Clocal_2D);
    Cv_in_2D = accumarray(rv, Clocal_2D);
    %Cv_2D = Cv_in_2D - Cv_out_2D;
    
    % get stationary distribution pi
    [W,~] = eigs(sto_mat',1);
    piDist = W/sum(W);
    
    % store results
    sampEntropy(iSamp).piDist = piDist;
    sampEntropy(iSamp).Slocal = Slocal;
    
    %sampCurvature_1D(iSamp).Clocal = Clocal_1D; % edge curvature
    sampCurvature_1D(iSamp).Cv_out = Cv_out_1D;
    sampCurvature_1D(iSamp).Cv_in = Cv_in_1D;
    
    sampCurvature_2D(iSamp).Cv_out = Cv_out_2D;
    sampCurvature_2D(iSamp).Cv_in = Cv_in_2D;
    
    toc(ti)
end
t2=toc(t1);
fprintf('\nDone with all %i samples\nTotal time elapsed: %0.2f sec\n',nSamp,t2)

% close parpool
delete(gcp('nocreate'))

%% save data
% compression?

fprintf('\nSaving data to %s\n',output_dir)
% make a new folder for this experiment
if ~exist(output_dir, 'dir')
   mkdir(output_dir)
end

sampName = matlab.lang.makeUniqueStrings(sampName);
sampName = matlab.lang.makeValidName(sampName);

fprintf('\n Saving data as .csv\n')

A = adj;
T = array2table(A,'RowNames',genes_max);
writetable(T, fullfile(output_dir,'adj.csv'), 'WriteRowNames', true);

A = expr_max;
T = array2table(A,'RowNames',genes_max,'VariableNames',sampName);
writetable(T, fullfile(output_dir,'expr.csv'), 'WriteRowNames', true);

A = [sampCurvature_1D(:).Cv_out];
T = array2table(A,'RowNames',genes_max,'VariableNames',sampName);
writetable(T, fullfile(output_dir,'1Dcurv_out.csv'), 'WriteRowNames', true);

A = [sampCurvature_1D(:).Cv_in];
T = array2table(A,'RowNames',genes_max,'VariableNames',sampName);
writetable(T, fullfile(output_dir,'1Dcurv_in.csv'), 'WriteRowNames', true);

A = [sampCurvature_2D(:).Cv_out];
T = array2table(A,'RowNames',genes_max,'VariableNames',sampName);
writetable(T, fullfile(output_dir,'2Dcurv_out.csv'), 'WriteRowNames', true);

A = [sampCurvature_2D(:).Cv_in];
T = array2table(A,'RowNames',genes_max,'VariableNames',sampName);
writetable(T, fullfile(output_dir,'2Dcurv_in.csv'), 'WriteRowNames', true);

% also save piDist, can use to recompute global averages
A = [sampEntropy(:).piDist];
T = array2table(A,'RowNames',genes_max,'VariableNames',sampName);
writetable(T, fullfile(output_dir,'piDist.csv'), 'WriteRowNames', true);

% save local entropy as well
A = [sampEntropy(:).Slocal];
T = array2table(A,'RowNames',genes_max,'VariableNames',sampName);
writetable(T, fullfile(output_dir,'nodeEntropy.csv'), 'WriteRowNames', true);

fprintf('data saved. script complete.\n')

end