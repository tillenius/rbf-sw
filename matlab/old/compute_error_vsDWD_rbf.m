%%%%% SPLINE INTERPOLATION TO RBF GRID FROM LAT-LON
% hmod=load('dwd/hmodc6_0000_E000605.dat'); % day zero
% hmod=load('dwd/hmodc6_0120_E000605.dat'); % day five
% hmod=load('dwd/hmodc6_0240_E000605.dat'); % day 10
% hmod=load('dwd/hmodc6_0360_E000605.dat'); % day 15

lad = hmod(:,1); thd = hmod(:,2); hd = hmod(:,3);
% lad(lad>180) = lad(lad>180)-360;
lad = lad*pi/180; thd = thd*pi/180;

la = atm.pts.la; th = atm.pts.th;
las = la; las(las<0)=las(las<0)+2*pi;
he = (H(:,4)+atm.gh0)/atm.g;
wghts = atm.pts.wghts;

% hsem=load('sem/williamson_6_0120_n16r5.dat');
% x = hsem(:,1); y = hsem(:,2); z = hsem(:,3);
% [la,th] = cart2sph(x,y,z);
% las = la; las(las<0) = las(las<0)+2*pi;
% he = hsem(:,4);
% wghts = hsem(:,end)';

% load('swtc6\p240.mat');
% la = p240(:,1); th = p240(:,2); he = p240(:,3);
% load('swtc6\gl240.mat');
% wghts = gl240(:,3).';
% load('swtc6\p120.mat');
% la = p120(:,1); th = p120(:,2); he = p120(:,3);
% load('swtc6\gl120.mat');
% wghts = gl120(:,3).';
% % la(la>180) = la(la>180)-360;
% la = la*pi/180; th = th*pi/180;
% las = la;

% h = latlon_interp2(lad,thd,hd,la,th);
h = latlon_interp2(lad,thd,hd,las,th);

diff = he-h;
% tri = delaunay(la,th);

errs = [wghts*abs(diff) sqrt(wghts*diff.^2) norm(diff,inf)]./...
    [wghts*abs(h) sqrt(wghts*h.^2) norm(h,inf)];

% xe = atm.pts.x; ye = atm.pts.y; ze = atm.pts.z;
% tri_l = TriangulateSpherePoints([xe ye ze]);
% wghts_l = ComputeTriQuadWeights(tri_l,[xe ye ze]);
% trisa_l = ComputeTriSurfaceArea(tri_l,[xe ye ze]);
% errs = [sum(trisa_l.*sum(wghts_l.*abs(diff(tri_l)),2)) ...
%         sqrt(sum(trisa_l.*sum(wghts_l.*diff(tri_l).^2,2))) ...
%         norm(diff,inf)]./...
%         [sum(trisa_l.*sum(wghts_l.*abs(he(tri_l)),2)) ...
%         sqrt(sum(trisa_l.*sum(wghts_l.*he(tri_l).^2,2))) ...
%         norm(he,inf)];