% PPI Network Topological Analysis
% By Kevin Murgas
% Stony Brook University, Dept. Biomedical Informatics

% loads in pre-processed PPI topology (vertices and edges)

% 1D simple graph model:
% compute # vertex, edges, and Euler characteristic
% Forman-Ricci curvature of only topology, using hop-metric (unit weights)
% other topological features can be examined as well (diameter, centrality)

% 2D simplicial complex model:
% define faces as triplets of vertices in feedforward or feedback direction
% compute # vertices, edges, and faces, Euler characteristic, FR curvature


%% load in data
% topo_stringdb11sparse
% derived from STRINGdb v11, with sparsification based on GO compartment localization terms
% relevant paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2652875/
datadir = 'Topology';
filename = 'topo_stringdb11sparse.mat';
load(fullfile(datadir,filename));
fprintf("Loaded network file: %s\n",filename)


%% 1D model: topological analysis
% print # vertices/edges, Euler characteristic

adj = full(adj);
nv = size(adj,1);
ne = sum(adj(:));
if issymmetric(adj) % if network is undirected (symmetric), then divide edges by 2
    ne = ne/2;
end

% euler characteristic: X = V - E + F
euler = nv - ne + 0;

fprintf(" 1-D topological features:\n #Vertex= %i, #Edge= %i, EulerChar= %i\n", nv, ne, euler)


%% 2D model: pre-analysis
% define triangular faces

% approach: undirected graph
% take subgraph of vertex neighbors (without vertex i)
% any edges present indicate a triangle, can store as a triple
% repeat but ignore previous vertices (already considered)

fprintf("Triangle finding: undirected graph method...")
adju = triu(adj); % take upper triangle of adj as undirected graph
[lv, rv] = ind2sub(size(adj), find(adju)); % vertex indices, only considering upper triangle (undirected)
nbrs = accumarray(lv, rv,[nv 1], @(x) {x}); % accumulate cell array of outgoing neighbors at each vertex

f_cell = cell(nv,1); % store faces as triplets in cell array
tic
for i=1:nv % scan through each vertex
    ni = nbrs{i}; % get adjacent neighbors (pre-computed)
    subadj = adju(ni,ni); % take subgraph of neighbors

    if any(subadj,'all') % any edges in subgraph indicate triangle
        [left, right] = ind2sub(size(subadj), find(subadj)); % find vertex pairs of edges
        faces = [left*0+i ni(left) ni(right)]; % create a list of triplets (3 columns: i,j,k)
        f_cell{i} = faces; % store in the cell of faces for this vertex
    else
        f_cell{i} = zeros(0,3); % otherwise store empty 0x3 array placeholder
    end
end
toc

% count total number of faces
nface_e = cellfun('size', f_cell, 1);
nface = sum(nface_e);


%% 2D model: topological analysis
% print # vertices/edges/faces, Euler characteristic

% euler characteristic: X = V - E + F
euler2 = nv - ne + nface;
fprintf("\n 2-D topological features:\n #Vertex= %i, #Edge= %i, #Face= %i, EulerChar= %i\n", nv, ne, nface, euler2)


%% Figure 3: edge and face degree distributions
% number of edges or faces incident to each vertex

% edge degree distribution:
deg = sum(adj,2); % out-going edge degree

% face degree distribution:
% count number of faces containing each vertex
all_face = cell2mat(f_cell);
counts = tabulate(all_face(:));
face_degree = 0*(1:nv);
face_degree(counts(:,1)) = counts(:,2);


figure(3);clf

% edge degree
subplot(1,2,1)
[N,edges] = histcounts(deg); % use histogram to bin before log-log transform
x = edges(2:end) - (edges(2)-edges(1))/2;
y = N/sum(N);
plot(x,y, 'ko','MarkerFaceColor','b','MarkerEdgeColor','none');
set(gca,'XScale','log','YScale','log') % set log-log scale
set(gca,'box','off','TickDir','out')
xlabel("edge.degree")
ylabel("P")
%title("Edge Degree Distribution")
axis square

% add log-log fit
hold on
xl=log(x); yl=log(y);
ind = y>0; % look where bin counts are non zero
[p,S] = polyfit(xl(ind),yl(ind),1); % linear fit with log(N) = b*log(k) + a  corresponding to power law fit N = exp(b * log(k) + a) = k^b * exp(a)
b = p(1);
a = exp(p(2));
R2 = 1 - (S.normr/norm(yl(ind) - mean(yl(ind))))^2;
yfit = a*x.^b;
plot(x,yfit, 'r')
hold off
text(x(10),yfit(3),sprintf("slope = %0.3f\nr^2 = %0.3f",b,R2))
xlim([5 1000])
ylim([10^-4.2 1])

%face degree
subplot(1,2,2)
[N,edges] = histcounts(face_degree);
x = edges(2:end) - (edges(2)-edges(1))/2;
y = N/sum(N);
plot(x,y, 'ko','MarkerFaceColor','b','MarkerEdgeColor','none');
set(gca,'XScale','log','YScale','log') % set log-log scale
set(gca,'box','off','TickDir','out')
xlabel("face.degree")
ylabel("P")
%title("Face Degree Distribution")
axis square

% add log-log fit
hold on
xl=log(x); yl=log(y);
ind = y>0; % look where bin counts are non zero
[p,S] = polyfit(xl(ind),yl(ind),1); % linear fit with log(N) = b*log(k) + a  corresponding to power law fit N = exp(b * log(k) + a) = k^b * exp(a)
b = p(1);
a = exp(p(2));
R2 = 1 - (S.normr/norm(yl(ind) - mean(yl(ind))))^2;
yfit = a*x.^b;
plot(x,yfit, 'r')
hold off
text(x(10),yfit(3),sprintf("slope = %0.3f\nr^2 = %0.3f",b,R2))
xlim([500 100000])
ylim([10^-4.2 1])


%% count triangles and quadrangles with matrix math
% exponentiate adjacency matrix for counting cycles/paths
A = sparse(adj)>0;
A2 = A^2;
A3 = A^3;
A4 = A2^2;

% get num of 2-cycles, 3-cycles, 4-cycles from diagonals
n2cyc = sum(diag(A2));
n3cyc = sum(diag(A3));
n4cyc = sum(diag(A4));

% number of unique triangles is 3-cycles divided by 6
% 3 start points * 2 directions = 6 path counts (repeated)
fprintf("# triangles = %i\n", full(n3cyc)/6);

% correct 4-cycles is only unique sets of 4 vertices: a->b->c->d->a
% remove pairs of 2-cycles (abaca) and 2 step loops (abcba)
n2cyc_pairs = sum(A2(:));
A2_nocyc = A2 - diag(diag(A2));
n2step = sum(A2_nocyc(:));
n4cyc_correct = n4cyc - n2cyc_pairs - n2step;

% number of unique quadangles is 4-cycles divided by 8
% 4 start points * 2 directions = 8 path counts (repeated)
fprintf("# quadrangles = %i\n", full(n4cyc_correct)/8);

%% count tetrahedrons
% 4 vertices that all share edges
% loop through all vertices, count if they all have edges

% pre-compute neighbors at each vertex
[lv, rv] = find(A);
nbrs = accumarray(lv, rv,[nv 1], @(x) {x'});

% begin count, storing hashed index
nct=0;
tet_hash = [];
b = nv+1; % base for index hash

tic
for i=1:nv
    fprintf("Vertex %i/%i (nct=%i)\n",i,nv,nct)
    for j=nbrs{i}(nbrs{i}>i)
        for k=nbrs{j}(nbrs{j}>j)
            if ~A(i,k)
                continue
            end
            for l=nbrs{k}(nbrs{k}>k)
                if A(i,l) && A(j,l)
                    nct = nct + 1;
                    %tet_hash = [tet_hash i*b^3 + j*b^2 + k*b + l];
                end
            end
        end
    end
end
toc

