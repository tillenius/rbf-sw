function mt_plot3d(par, H)

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

DT = delaunayTriangulation(nodes(:,1), nodes(:,2), nodes(:,3));
hullFacets = convexHull(DT);


clf;

colormap('jet')

axis equal

%%% set viewpoint
global az el
az=50;
el=20;

% rotate object instead of changing viewpoint to avoid resize
%view(az, el)
nodes = nodes*rotz(az)*rotx(-el);
view(0,0)

axis off
hold on

%%% draw circle around sphere

% find view direction
%[az,el]=view
%n = rotz(az)*rotx(-el)*[0 1 0]';
%circ3d(n)
circ3d([0 1 0]');

%%% draw solution

hold on

trisurf(hullFacets, nodes(:,1), nodes(:,2), nodes(:,3), zeta, 'EdgeColor', 'none');

switch par.test
  case 'tc5' % Draw the mountain

    % Parameters for the mountain:
    % cone of radius mR=pi/9 located at (lam_c, thm_c), height hm0=2000
    atm.lam_c = -pi/2;
    atm.thm_c = pi/6;
    atm.mR = pi/9;
    atm.hm0 = 2000;

    theta=0:0.01:2*pi;
    plotsph(cos(theta)*atm.mR + atm.lam_c, sin(theta)*atm.mR + atm.thm_c, 'w:')
end


%%% draw some longitudes
theta=0:0.01:2*pi;

plotsph( theta,  0/90*pi/2*ones(size(theta)), 'k-')
for lat=[-67.5 -45 -22.5 22.5 45 67.5 90]
  plotsph( theta,  lat/90*pi/2*ones(size(theta)), 'k--')
end

plotsph( 0/90*pi/2*ones(size(theta)), theta, 'k-')
for long=[-67.5 -45 -22.5 22.5 45 67.5 90]
  plotsph( long/90*pi/2*ones(size(theta)), theta, 'k--')
end


%%%%%%%%%%%%%%%

function Rx = rotx(a)
  Rx = [ 1 0 0 ; 0 cosd(a) -sind(a) ; 0 sind(a) cosd(a) ];

function Ry = roty(a)
  Ry = [ cosd(a) 0 sind(a) ; 0 1 0 ; -sind(a) 0 cosd(a) ];

function Rz = rotz(a)
  Rz = [ cosd(a) -sind(a) 0 ; sind(a) cosd(a) 0 ; 0 0 1 ];

function circ3d(n)
  theta=0:0.01:2*pi;
  v=null(n');
  points=repmat([0 0 0]',1,size(theta,2))+1*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
  plot3(points(1,:),points(2,:),points(3,:),'k-','linewidth',4);

function plotsph(lam, thm, sty)
global az el
%    pts = [cos(theta)*atm.mR + atm.lam_c ; sin(theta)*atm.mR + atm.thm_c];
%    [x,y,z] = sph2cart(pts(1,:), pts(2,:), 1);
    [x,y,z] = sph2cart(lam, thm, 1);
    A=[x',y',z']*rotz(az)*rotx(-el);
    plot3(A(:,1),A(:,2),A(:,3),sty,'linewidth',2);
