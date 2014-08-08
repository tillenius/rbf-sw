function [ind_i, ind_j, weightsDx, weightsDy, weightsDz, weightsL] = mt_rbfmatrix_fd_hyper(x,tree,par,atma)
%%% [D,L] = rbfmatrix_fd_tree(x,ep,alpha,fdsize,order,dim)
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

rbf = @(ep,rd2) exp(-ep^2*rd2);
drbf = @(ep,rd2) -2*ep^2*exp(-ep^2*rd2);

weightsDx = zeros(N*fdsize,1);
weightsDy = zeros(N*fdsize,1);
weightsDz = zeros(N*fdsize,1);
weightsL = zeros(N*fdsize,1);

ind_i = zeros(N*fdsize,1);
ind_j = zeros(N*fdsize,1);

%A = ones(fdsize+1,fdsize+1); A(end,end) = 0;
%B = zeros(fdsize+1,4);

global cnd
global sz
cnd = zeros(N,1);
sz = zeros(N,1);

tree = knnsearch(x,x,'k',fdsize)';

idxctr = 1;
c = par.gamma*ep^(2*order);
disp('loop')
for j=1:N

    badness = 1e12;
    cfdsize = fdsize;
    first = true;
    while badness > 3e11
        imat = tree(:, j);
        imat = imat(1:cfdsize,:);
        
        rd2 = max(0,2*(1-x(imat,1)*x(imat,1).'-...
            x(imat,2)*x(imat,2).'-x(imat,3)*x(imat,3).'));
        
        A = ones(cfdsize+1,cfdsize+1); A(end,end) = 0;
        A(1:cfdsize,1:cfdsize) = rbf(ep,rd2);

        badness = cond(A);

%        if ~first
%            disp([num2str(j) ': cfdsize=' num2str(cfdsize) ' badness=' num2str(badness)]);
%        end

        cfdsize = cfdsize - 1;
        first = false;
    end
    cfdsize = cfdsize + 1;

    outidx=idxctr:(idxctr+cfdsize-1);
    dp = (x(j,1)*x(imat,1) + x(j,2)*x(imat,2) + x(j,3)*x(imat,3));
    rd2v = rd2(:,1);
    dr = drbf(ep,rd2v);
    B = zeros(cfdsize+1,4);
    B(1:cfdsize,1:3) = [(x(j,1)*dp-x(imat,1)).*dr, ...
        (x(j,2)*dp-x(imat,2)).*dr, (x(j,3)*dp-x(imat,3)).*dr]/atma;
    B(1:cfdsize,4) = c*hyper(ep^2*rd2v,dim,order).*exp(-ep^2*rd2v);
    weights = A\B;

    cnd(j) = badness;
    sz(j) = cfdsize;

    t = sortrows([imat weights(1:cfdsize,:)], 1);

    ind_i(outidx) = j;
    ind_j(outidx) = t(:,1);

    weightsDx(outidx) = t(:,2);
    weightsDy(outidx) = t(:,3);
    weightsDz(outidx) = t(:,4);
    weightsL(outidx) = t(:,5);
    idxctr = idxctr + cfdsize;
end

ind_i = ind_i(1:idxctr-1);
ind_j = ind_j(1:idxctr-1);

weightsDx = weightsDx(1:idxctr-1);
weightsDy = weightsDy(1:idxctr-1);
weightsDz = weightsDz(1:idxctr-1);
weightsL  = weightsL(1:idxctr-1);

function p=hyper(ep2r2,d,k)
%%% laplacian to the power k of dimension d
% ep2r2 - ep^2*r^2
n = length(ep2r2);
P = zeros(n,k+1);
P(:,1) = 1; P(:,2) = 4*ep2r2-2*d;
for j=3:k+1
    P(:,j) = 4*(ep2r2-2*j-d/2+4).*P(:,j-1) - ...
            8*(j-2)*(2*j+d-6)*P(:,j-2);
end
p = P(:,k+1);
