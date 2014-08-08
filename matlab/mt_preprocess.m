function mt_preprocess(nodes, fdsize)

if ischar(nodes)
    nodes=load(nodes);
    if isstruct(nodes)
        nodes=nodes.nodes;
    end
end

[la, th] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
[x,y,z] = sph2cart(la, th, ones(size(la)));
nodes = [x,y,z];
clear x y z

[nodes, neighbors] = mt_permute(nodes, fdsize);

numnodes = length(nodes);

nodesname = ['../nodesets/' num2str(numnodes) '_' num2str(fdsize) '.mat'];
treename  = ['../tree/'     num2str(numnodes) '_' num2str(fdsize) '.mat'];
disp(['save nodes "' nodesname '"']);
save(nodesname, 'nodes');
disp(['save tree "' treename '"']);
save(treename, 'neighbors');



function [nodes, neighbors] = mt_permute(nodes, fdsize)

N = length(nodes);
ind_i = zeros(N*fdsize,1);
ind_j = zeros(N*fdsize,1);

neighbors = knnsearch(nodes,nodes,'k',fdsize)';

for j=1:N
    imat = neighbors(:,j);
    ind_i((j-1)*fdsize+1:j*fdsize) = j;
    ind_j((j-1)*fdsize+1:j*fdsize) = imat;
end

A = sparse(ind_i, ind_j, 1);
perm = symrcm(A);

nodes = nodes(perm,:);

% update indices according to permutation
iperm(perm) = 1:N;
neighbors = iperm(neighbors);

% update column order according to permutation
neighbors = neighbors(:,perm);
