function [ind_i,ind_j,weightsDx,weightsDy,weightsDz] = mt_rbfmatrix_fd_cart(x,tree,fdsize,ep)
%%% [Dx,Dy,Dz] = rbfmatrix_fd_cart(x,ep,alpha,fdsize)
% Requires kd-tree code.
% IN:
% x - nodes
% ep - shape parameter
% fdsize - stencil size (good choices: 31, 50, 74, 101)
% OUT:
% Dx - sparse differentiation matrix
% Dy - sparse differentiation matrix
% Dy - sparse differentiation matrix

N = length(x);
% srange = sqrt(6*fdsize/N); % search range

rbf = @(ep,rd2) exp(-ep^2*rd2);
drbf = @(ep,rd2) -2*ep^2*exp(-ep^2*rd2);

weightsDx = zeros(N*fdsize,1);
weightsDy = zeros(N*fdsize,1);
weightsDz = zeros(N*fdsize,1);

ind_i = zeros(N*fdsize,1);
ind_j = zeros(N*fdsize,1);

A = ones(fdsize+1,fdsize+1); A(end,end) = 0;
B = zeros(fdsize+1,1);

disp('loop')
for j=1:N

    imat = tree(:, j);
    ind_i((j-1)*fdsize+1:j*fdsize) = j;
    ind_j((j-1)*fdsize+1:j*fdsize) = imat;
    
    dp = (x(j,1)*x(imat,1) + x(j,2)*x(imat,2) + x(j,3)*x(imat,3));
    
    rd2 = max(0,2*(1-x(imat,1)*x(imat,1).'-...
        x(imat,2)*x(imat,2).'-x(imat,3)*x(imat,3).'));
    rd2v = rd2(:,1);
    
    A(1:fdsize,1:fdsize) = rbf(ep,rd2);
    [LA,UA,P] = lu(A);

    B(1:fdsize) = (x(j,1) - x(imat,1)).*drbf(ep,rd2v);
    weights = UA\(LA\(P*B));
    weightsDx((j-1)*fdsize+1:j*fdsize) = weights(1:fdsize);

    B(1:fdsize) = (x(j,2) - x(imat,2)).*drbf(ep,rd2v);
    weights = UA\(LA\(P*B));
    weightsDy((j-1)*fdsize+1:j*fdsize) = weights(1:fdsize);

    B(1:fdsize) = (x(j,3) - x(imat,3)).*drbf(ep,rd2v);
    weights = UA\(LA\(P*B));
    weightsDz((j-1)*fdsize+1:j*fdsize) = weights(1:fdsize);

    
end
