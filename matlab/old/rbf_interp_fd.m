function he = rbf_interp_fd(x,ep,fd,h,xe)
% he = rbf_interp_fd(x,ep,fd,h,xe)
% Interpolate using RBF-FD
% This routine could use some improvement in terms of efficiency.

N = length(x);
Ne = length(xe);

% srange = sqrt(6*fd/N); % search range

rbf = @(ep,rd2) exp(-ep^2*rd2);

% [idx_e, dist, root] = kdtreeidx(x,xe);
he = zeros(Ne,1);

root = kdtree_build(x);
idx_e = kdtree_nearest_neighbor(root, xe);
[idx_e,idx_o] = sort(idx_e);
idx_current = 0;
%%% For every node in the evaluation node set, the closest node in the
%%% original node set is found and an interpolant from its stencil is used
%%% to evaluate at the new node location.
for j = 1:Ne
%     [pts,rdv,idx] = kdrangequery(root, x(idx_e(j),:), srange);
%     [tmp,ir] = sort(rdv);
%     ir = ir(1:fd);
 
%     imat = idx(ir);
    if idx_current ~= idx_e(j)
        idx = kdtree_k_nearest_neighbors(root, x(idx_e(j),:), fd);
        imat = idx(fd:-1:1);
        xi = x(imat,:);
        
        rd2 = max(0,2*(1-xi(:,1)*xi(:,1).'-...
            xi(:,2)*xi(:,2).'-xi(:,3)*xi(:,3).'));
        
        A = rbf(ep,rd2);
%         [LA,UA,P] = lu(A);
%         lambda = UA\(LA\(P*hi));
        lambda = A\h(imat);
    end
%     rd2e = max(0,2*(1-xe(j,1)*x(imat,1).'-...
%         xe(j,2)*x(imat,2).'-xe(j,3)*x(imat,3).'));
%     rd2e = max(0,2*(1-xe(idx_o(j),1)*x(imat,1).'-...
%         xe(idx_o(j),2)*x(imat,2).'-xe(idx_o(j),3)*x(imat,3).'));
    rd2e = max(0,2*(1-xe(idx_o(j),1)*xi(:,1)-...
        xe(idx_o(j),2)*xi(:,2)-xe(idx_o(j),3)*xi(:,3))).';
    Ae = rbf(ep,rd2e);
    he(idx_o(j)) = Ae*lambda;
    idx_current = idx_e(j);
end
% kdtree([], [], root);
kdtree_delete(root);