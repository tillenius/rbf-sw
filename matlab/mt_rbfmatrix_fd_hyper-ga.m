function [ind_i, ind_j, weightsDx, weightsDy, weightsDz, weightsL] = rbfmatrix_fd_hyper_proj_ga_choose(x,tree,par,atma)
%%% [D,L] = rbfmatrix_fd_tree(x,ep,alpha,fdsize,order,dim)
% Requires kd-tree code.
% IN:
% x - nodes
% ep - shape parameter
% fdsize - stencil size (good choices: 31, 50, 74, 101)
% order - L = L^order
% dim - dimension of Laplacian formula
% OUT:
% DPx - sparse differentiation matrix
% DPy - sparse differentiation matrix
% DPy - sparse differentiation matrix
% L - sparse dissipation matrix

N = length(x);
ep = par.ep;
fdsize = par.fd;
order = par.order;
dim = par.dim;

weightsDx = zeros(N*fdsize,1);
weightsDy = zeros(N*fdsize,1);
weightsDz = zeros(N*fdsize,1);
weightsL = zeros(N*fdsize,1);

ind_i = zeros(N*fdsize,1);
ind_j = zeros(N*fdsize,1);

warning('OFF', 'MATLAB:nearlySingularMatrix')
msg_length = 0;
for j=1:N

    imat = tree(:, j);
    outidx=(j-1)*fdsize+1:j*fdsize;

    xi = x(imat,:);
    P = [(x(j,2)^2+x(j,3)^2), -x(j,1)*x(j,2), -x(j,1)*x(j,3);
        -x(j,1)*x(j,2), (x(j,1)^2+x(j,3)^2), -x(j,2)*x(j,3);
        -x(j,1)*x(j,3), -x(j,2)*x(j,3), (x(j,1)^2+x(j,2)^2)];
    xp = repmat(x(j,:),fdsize,1)+xi*P;
    
    %weights = zeros(fdsize,4);


    opts.const = 1;
    opts.verbose = 0;
    in = 1:fdsize;
    wxyz = rbfga_weights({'x','y','z'},ep,xp(in,:),xp(1,:),opts);
    wL = par.gamma*rbfga_weights(['L' num2str(order)],ep,xp(in,:),xp(1,:),opts);

    t = sortrows([imat wxyz(1:fdsize,:) wL(1:fdsize,:)], 1);

    ind_i(outidx) = j;
    ind_j(outidx) = t(:,1);

    weightsDx(outidx) = t(:,2)/atma;
    weightsDy(outidx) = t(:,3)/atma;
    weightsDz(outidx) = t(:,4)/atma;
    weightsL(outidx) = t(:,5);
    
end
%fprintf(1,'\n');
%DPx = sparse(I,J,weightsDx,N,N); 
%DPy = sparse(I,J,weightsDy,N,N); 
%DPz = sparse(I,J,weightsDz,N,N); 
%L = sparse(I,J,weightsL,N,N);
warning('ON', 'MATLAB:nearlySingularMatrix')
