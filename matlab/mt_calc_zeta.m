function zeta=mt_calc_zeta(par, H)

name=[ num2str(par.n) '_' num2str(par.fd) '.mat'];
epname=[ num2str(par.n) '_' num2str(par.fd) '_ep' num2str(par.ep) '.mat'];

nodes = getnodes(par);

try
  d = load(['../cartfd/' epname]);
catch
  disp(['MAKECART ']);
  par
  mt_makecart(par);
  d = load(['../cartfd/' epname]);
end

i=d.i; j=d.j; dx=d.dx; dy=d.dy; dz=d.dz;

X = H(:,1); Y = H(:,2); Z = H(:,3);
clear H

N = length(X);

Dx = sparse(i,j,dx,N,N); 
DxZ = Dx*Z;
DxY = Dx*Y;
clear Dx

Dy = sparse(i,j,dy,N,N); 
DyZ = Dy*Z;
DyX = Dy*X;
clear Dy

Dz = sparse(i,j,dz,N,N);
DzY = Dz*Y;
DzX = Dz*X;
clear Dz

atm.a = 6.37122e6;  % galew_setup.m

zeta = (nodes(:,1).*(DyZ-DzY) + nodes(:,2).*(DzX-DxZ) + nodes(:,3).*(DxY-DyX))/atm.a;
