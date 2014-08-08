function [C,h] = plotlatlong_fd(x,ep,fd,h,res,minh,maxh,levels)
% plotlatlong_fd(x,ep,fd,h,res,minh,maxh,levels)
la = linspace(-pi,pi,2*res);
th = linspace(-pi/2,pi/2,res);
% th = linspace(0,pi/2,res);
[L,T] = meshgrid(la,th);
[xe,ye,ze] = sph2cart(L(:),T(:),ones(2*res^2,1));
% ne = length(xe);
xe = [xe ye ze];

he = rbf_interp_fd(x,ep,fd,h,xe);
disp(['max: ' num2str(max(he),'%.2e') ', min: ' num2str(min(he),'%.2e') '.']);

H = reshape(he,res,2*res);
figure;
% [C,h] = contourf(L*360/(2*pi),T*360/(2*pi),H,linspace(minh,maxh,levels));
% [C,h] = contourf(L*360/(2*pi),T*360/(2*pi),H,...
%     [minh:levels:-levels levels:levels:maxh]);
% [C,h] = contour(L*360/(2*pi),T*360/(2*pi),H,...
%     [minh:levels:-levels levels:levels:maxh]);
% [C,h] = contourf(L*360/(2*pi),T*360/(2*pi),H,[minh:levels:maxh]);
% [C,h] = contour(L*360/(2*pi),T*360/(2*pi),H,[minh:levels:maxh]);

% surf(L*180/pi,T*180/pi,H); shading interp; view([0 90]); C=0; h=0;
% surf(L*180/pi,T*180/pi,zeros(size(L)),H); shading interp; view([0 90]); C=0; h=0;
surf(L*180/pi,T*180/pi,H,0.1*ones(size(L))); shading interp; C=0; h=0;
colormap(gray);
xlabel('Longitude','FontSize',12);
ylabel('Latitude','FontSize',12);
zlabel('Height','FontSize',12);
view([-24 88]); set(gca,'FontSize',11);
axis([-180 180 -90 90 0 2000]);
% xlabel('Longitude','Interpreter','latex','FontSize',12);
% ylabel('Latitude','Interpreter','latex','FontSize',12);
% colorbar;
% set(gca,'FontSize',11);
% axis([-180 180 -90 90]);
% axis([-180 180 0 90]);
% axis equal