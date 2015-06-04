function [DT, hullFacets] = mt_plot3(par, H, DT, hullFacets)

nodes = getnodes(par);

atm.g = 9.80616;       % Gravitational constant (m/s^2).

switch par.test
  case 'galew'
    zeta = mt_calc_zeta(par, H);
  case 'tc5'
    atm.gh0 = atm.g*5960;     % Initial condition for the geopotential field (m^2/s^2).
    zeta = (H(:,4)+atm.gh0)/atm.g;
  otherwise
    error(['Unknown test case "' par.test '"']);
end

if nargin ~= 4
  DT = delaunayTriangulation(nodes(:,1), nodes(:,2), nodes(:,3));
  hullFacets = convexHull(DT);
end


clf;

colormap('jet')

axis equal
axis([-1 1 -1 1 -1 1])

% set viewpoint
global az el 
az=50;
el=20;

% rotate object instead of changing viewpoint to avoid resize
nodes = nodes*rotz(az)*rotx(-el);
view(0,0)

axis off
hold on

%%% draw solution

H = trisurf(hullFacets, nodes(:,1), nodes(:,2), nodes(:,3), zeta, 'EdgeColor', 'none');
shading interp;

%%% draw a black circle disc through the sphere to get an outline around it

circ3d([0 1 0]');


%%% draw some latitudes and longitudes

theta=0:0.01:2*pi;

%lats=[-67.5 -45 -22.5 22.5 45 67.5 90];
lats=[-60 -30 30 60 90];
plotsph( theta,  0/90*pi/2*ones(size(theta)), 'k-', 1)
for lat=lats
  plotsph( theta,  lat/90*pi/2*ones(size(theta)), 'k--', 1)
end

plotsph( 0/90*pi/2*ones(size(theta)), theta, 'k-', 1)
for long=lats
  plotsph( long/90*pi/2*ones(size(theta)), theta, 'k--', 1)
end


%%% add a mountain for the tc5 case

switch par.test
  case 'tc5' % Draw the mountain

    % Parameters for the mountain:
    % cone of radius mR=pi/9 located at (lam_c, thm_c), height hm0=2000
    atm.lam_c = -pi/2;
    atm.thm_c = pi/6;
    atm.mR = pi/9;
    atm.hm0 = 2000;

    theta=0:0.005:2*pi;
    plotsph(cos(theta)*atm.mR + atm.lam_c, sin(theta)*atm.mR + atm.thm_c, 'k-', 2)

end


%%% set window size, position and camera and take a screenshot.

set(gcf, 'Position', [400, 150, 800, 800]) % set window size to 800 x 800
camva(5.4);                        % zoom in to fill window
a = gca;
a.Position = [0.095 .095 .81 .81]; % move sphere to center and adjust zoom
set(gcf,'color','w');              % white background instead of gray
img = getframe(gcf);               % get rgb data
imwrite(img.cdata, 'bitmap.png');  % write rgb data to file


%%%%%%%%%%%%%%%

function Rx = rotx(a)
  Rx = [ 1 0 0 ; 0 cosd(a) -sind(a) ; 0 sind(a) cosd(a) ];

function Ry = roty(a)
  Ry = [ cosd(a) 0 sind(a) ; 0 1 0 ; -sind(a) 0 cosd(a) ];

function Rz = rotz(a)
  Rz = [ cosd(a) -sind(a) 0 ; sind(a) cosd(a) 0 ; 0 0 1 ];

function circ3d(n)
  global x y z
  theta=0:0.005:2*pi;
  v=null(n');
  points=repmat([0 0 0]',1,size(theta,2))+1*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
  patch(points(1,:), points(2,:), points(3,:), 'k');

function plotsph(lam, thm, sty, lsize)
  global az el
  [x,y,z] = sph2cart(lam, thm, 1);
  A=[x',y',z']*rotz(az)*rotx(-el);
  plot3(A(:,1),A(:,2),A(:,3), sty, 'LineWidth', lsize);
