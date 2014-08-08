% Initial condition for test case 5.  uc contains the velocity in Cartesian 
% coordinates with uc(:,1:3) the x,y,z direction, respectively.  gh contains
% the geopotential height field from the mean geopotential gh0.
function [uc,gh,us] = tc5_computeInitialCondition(atm,t)

% Extract out the constants from the atm structure to make the code easier
% to read.
alpha = atm.alpha; a = atm.a; omega = atm.omega; u0 = atm.u0;
x = atm.pts.x; y = atm.pts.y; z = atm.pts.z;

gh = - (a*omega*u0+ u0^2/2)*(-x*sin(alpha) + z*cos(alpha)).^2;
uc = u0*[-y*cos(alpha) (x*cos(alpha) + z*sin(alpha)) -y*sin(alpha)];

% vectors for translating the field in Cartesian coordinates to a field
% in spherical coordinates.
c2s_u = [-sin(atm.pts.la) -sin(atm.pts.th).*cos(atm.pts.la)];
c2s_v = [cos(atm.pts.la) -sin(atm.pts.th).*sin(atm.pts.la)];
c2s_w = [zeros(size(atm.pts.la)) cos(atm.pts.th)];

us(:,1) = c2s_u(:,1).*uc(:,1) + c2s_v(:,1).*uc(:,2) + c2s_w(:,1).*uc(:,3);
us(:,2) = c2s_u(:,2).*uc(:,1) + c2s_v(:,2).*uc(:,2) + c2s_w(:,2).*uc(:,3);
