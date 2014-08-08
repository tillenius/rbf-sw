function [Dx,Dy,Dz] = rbfmatrix_global_cart(xc,ep)
% Create differentiation matrices using global RBFs
% INPUT:
% xc - node locations in Cartesian coordinates
% ep - shape parameter
% OUTPUT:
% Dx, Dy, Dz - Cartesian differentiation matrices

rbf = inline('1./sqrt(1 + ep^2.*r)'); % IMQ
drbfor = inline('-ep^2./(sqrt(1 + ep^2.*r).^3)'); 
nd = length(xc);
[xk,xj] = meshgrid(xc(:,1));
[yk,yj] = meshgrid(xc(:,2));
[zk,zj] = meshgrid(xc(:,3));

dp = (xj.*xk + yj.*yk + zj.*zk);

rd2 = 2 - 2*dp;

B = drbfor(ep,rd2);

Dx = (xj - xk).*B;
Dy = (yj - yk).*B;
Dz = (zj - zk).*B;

clear B xj xk yj yk zj zk temp dp;

A = rbf(ep,rd2);

opts.SYM = true;
% % % opts.POSDEF = true;
disp('Calculating Dx...');
Dx = linsolve(A,Dx.',opts).';
disp('Calculating Dy...');
Dy = linsolve(A,Dy.',opts).';
disp('Calculating Dz...');
Dz = linsolve(A,Dz.',opts).';

% temp = [Dx;Dy;Dz]/A;
% Dx = temp(1:nd,:);
% Dy = temp(nd+1:2*nd,:);
% Dz = temp(2*nd+1:3*nd,:);
