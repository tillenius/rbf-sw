function plotlatlong(x,y,z,ep,h,res)
% plotlatlong(x,y,z,ep,h,res)
% Plot a filled contour plot on a lat/long grid
% This uses global RBF interpolation
la = linspace(-pi,pi,2*res);
th = linspace(-pi/2,pi/2,res);
[L,T] = meshgrid(la,th);
[xe,ye,ze] = sph2cart(L(:),T(:),ones(2*res^2,1));
ne = length(xe);

% rbf = @(ep,rd2) sqrt(1+ep^2*rd2);
rbf = @(ep,rd2) 1./sqrt(1+ep^2*rd2);
rd2 = max(2*(1-x*x.'-y*y.'-z*z.'),0);
A = rbf(ep,rd2);
coef = A\h;

he = zeros(ne,1);
for j=1:ne
    rd2 = max(2*(1-xe(j)*x.'-ye(j)*y.'-ze(j)*z.'),0);
    A = rbf(ep,rd2);
    he(j) = A*coef;
end
disp(['max: ' num2str(max(he),'%.2e') ', min: ' num2str(min(he),'%.2e') '.']);

H = reshape(he,res,2*res);
figure;
% contourf(L*360/(2*pi),T*360/(2*pi),H,linspace(5e3,6e3,21));
% contourf(L*360/(2*pi),T*360/(2*pi),H,linspace(8e3,11e3,16));
contourf(L*360/(2*pi),T*360/(2*pi),H,...
    [-12e-5:2e-5:-2e-5 2e-5:2e-5:16e-5]);

% surf(L*360/(2*pi),T*360/(2*pi),H), shading interp; view([0 90]);
axis([-180 180 -90 90]);
axis([-180 180 0 90]);
xlabel('Longitude','Interpreter','latex','FontSize',12);
ylabel('Latitude','Interpreter','latex','FontSize',12);
colorbar;
set(gca,'FontSize',11);