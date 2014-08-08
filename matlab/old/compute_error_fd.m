% %% Compare using RBF-FD interpolation

% load('data/h163842.mat') % l2 error 5e-4 compared to DWD
% xe = atme.pts.x; ye = atme.pts.y; ze = atme.pts.z;
% la = atme.pts.la; th = atme.pts.th;
% wghts = atme.pts.wghts;
% ne = length(xe);

% load('swtc6\p024.mat'); p = p024; % day one
% load('swtc6\p120.mat'); p = p120; % day five
% load('swtc6\p168.mat'); p = p168; % day seven
% load('swtc6\p336.mat'); p = p336; % day fourteen
% la = p(:,1); th = p(:,2); h = p(:,3);
% load('swtc6\h120.mat');
% load('swtc6\gl120.mat'); wghts = gl120(:,3)';
% load('swtc6\h240.mat');
% load('swtc6\gl240.mat'); wghts = gl240(:,3)';
% load('swtc6\h336.mat');
% load('swtc6\gl336.mat'); wghts = gl336(:,3)';
% la(la>180) = la(la>180)-360;
% la = la*pi/180; th = th*pi/180;
% ne = length(la);
% % [xe,ye,ze] = sph2cart(la*2*pi/360,th*2*pi/360,ones(ne,1));
% [xe,ye,ze] = sph2cart(la,th,ones(ne,1));

% hsem=load('sem/williamson_6_0120_n16r5.dat');
% xe = hsem(:,1); ye = hsem(:,2); ze = hsem(:,3);
% [la,th] = cart2sph(xe,ye,ze);
% h = hsem(:,4);
% wghts = hsem(:,end)';

% refname = 'data\TC5_n25600_ep9.12_fd101g';
refname = 'data/TC5_n163842_ep14.0671_fd31g'; 
load([refname '.mat'], 'H', 'atm');
xe = atm.pts.x; ye = atm.pts.y; ze = atm.pts.z;
la = atm.pts.la; th = atm.pts.th;
h = (H(:,4)+atm.gh0)/atm.g;
wghts = atm.pts.wghts;
load([opts.bsname '.mat']);

tri = delaunay(la,th);
x = atm.pts.x; y = atm.pts.y; z = atm.pts.z;

he = rbf_interp_fd([x y z],par.ep,par.fd,(H(:,4)+atm.gh0)/atm.g,[xe ye ze]);

diff = he-h;
errs = [wghts*abs(diff) sqrt(wghts*diff.^2) norm(diff,inf)]./...
    [wghts*abs(h) sqrt(wghts*(h).^2) norm(h,inf)];
% errs = [norm(diff,1) norm(diff) norm(diff,inf)]./...
%     [norm(h,1) norm(h,2) norm(h,inf)]; % No quadrature weights
% errs = [sum(trisa.*sum(wghts.*abs(diff(trie)),2)) ...
%     sqrt(sum(trisa.*sum(wghts.*diff(trie).^2,2))) norm(diff,inf)]./...
%     [sum(trisa.*sum(wghts.*abs(h(trie)),2)) ...
%     sqrt(sum(trisa.*sum(wghts.*h(trie).^2,2))) norm(h,inf)];