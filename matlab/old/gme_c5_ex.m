% Load in the data

% load('swtc6\p168.mat'); % day seven
% la = p168(:,1); th = p168(:,2); h = p168(:,3);
load('swtc6\p336.mat'); % day fourteen
la = p336(:,1); th = p336(:,2); h = p336(:,3);
la(la>180) = la(la>180)-360;
ne = length(la);
[xe,ye,ze] = sph2cart(la*2*pi/360,th*2*pi/360,ones(ne,1));

tri = delaunay(la,th);
x = atm.pts.x; y = atm.pts.y; z = atm.pts.z;

rbf = @(ep,rd2) sqrt(1+ep^2*rd2);
rd2 = max(2*(1-x*x.'-y*y.'-z*z.'),0);
A = rbf(ep,rd2);
coef = A\((H(:,4)+atm.gh0)/atm.g);

he = zeros(ne,1);
for j=1:ne
    rd2 = max(2*(1-xe(j)*x.'-ye(j)*y.'-ze(j)*z.'),0);
    A = rbf(ep,rd2);
    he(j) = A*coef;
end

diff = he-h;
% errs = [wghts*abs(diff) sqrt(wghts*diff.^2) norm(diff,inf)]./...
%     [wghts*abs(h) sqrt(wghts*(h).^2) norm(h,inf)];
errs = [norm(diff,1) norm(diff) norm(diff,inf)]./...
    [norm(h,1) norm(h,2) norm(h,inf)]; % No quadrature weights