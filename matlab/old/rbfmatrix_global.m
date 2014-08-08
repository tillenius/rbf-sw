function [DPx,DPy,DPz,L] = rbfmatrix_global(xc,ep)
% Create differentiation matrices using global RBFs
% INPUT:
% xc - node locations in Cartesian coordinates
% ep - shape parameter
% OUTPUT:
% DPx, DPy, DPz - projected differentiation matrices
% L - A^-1.
rbf = inline('1./sqrt(1 + ep^2.*r)'); % IMQ
drbfor = inline('-ep^2./(sqrt(1 + ep^2.*r).^3)'); 
nd = length(xc);
[xk,xj] = meshgrid(xc(:,1));
[yk,yj] = meshgrid(xc(:,2));
[zk,zj] = meshgrid(xc(:,3));

dp = (xj.*xk + yj.*yk + zj.*zk);

rd2 = 2 - 2*dp;

B = drbfor(ep,rd2);

DPx = (xj.*dp - xk).*B;
DPy = (yj.*dp - yk).*B;
DPz = (zj.*dp - zk).*B;

clear B xj xk yj yk zj zk temp dp;

A = rbf(ep,rd2);

opts.SYM = true;
% % opts.POSDEF = true;
disp('Calculating DPx...');
DPx = linsolve(A,DPx.',opts).';
disp('Calculating DPy...');
DPy = linsolve(A,DPy.',opts).';
disp('Calculating DPz...');
DPz = linsolve(A,DPz.',opts).';
disp('Calculating L...');
L = linsolve(A,eye(nd),opts).';
disp('...done.');
% temp = [DPx;DPy;DPz;eye(nd)]/A;
% DPx = temp(1:nd,:);
% DPy = temp(nd+1:2*nd,:);
% DPz = temp(2*nd+1:3*nd,:);
% L = temp(3*nd+1:4*nd,:);