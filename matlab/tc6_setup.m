%
% Set up for the Williamson test case 6.
%
function atm = tc6_setup(nfile)

if ischar(nfile)
   nodes = load(nfile);
else
   nodes = nfile;
end
[atm.pts.la,atm.pts.th,r] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
[nodes(:,1),nodes(:,2),nodes(:,3)] = sph2cart(atm.pts.la,atm.pts.th,ones(size(atm.pts.la)));

atm.pts.nodes = nodes(:,1:3);
atm.pts.x = nodes(:,1); atm.pts.y = nodes(:,2); atm.pts.z = nodes(:,3);
atm.pts.nd = size(atm.pts.x,1);
if size(nodes,2) == 4   % Extract the quadrature weights.
   atm.pts.wghts = nodes(:,4)';
else
   atm.pts.wghts = repmat(atm.pts.nd/(4*pi),[1 atm.pts.nd]);
end

[atm.pts.la,atm.pts.th,r] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3)); 

% Variables for projecting an arbitrary Cartesian vector onto the surface
% of the sphere.
x2 = nodes(:,1).^2; xy = nodes(:,1).*nodes(:,2);
y2 = nodes(:,2).^2; xz = nodes(:,1).*nodes(:,3);
z2 = nodes(:,3).^2; yz = nodes(:,2).*nodes(:,3);
atm.pts.p_u = [1-x2  -xy   -xz];
atm.pts.p_v = [-xy  1-y2   -yz];
atm.pts.p_w = [-xz   -yz  1-z2];

atm.a = 6.37122e6;     % Mean radius of the earth (meters).
atm.alpha = 0;          % Angle of rotation
atm.omega = 7.848e-6;  % A constant (1/s)
atm.Omega = 7.292e-5; % Rotation rate of the earth (1/seconds).
atm.g = 9.80616;       % Gravitational constant (m/s^2).
atm.gh0 = atm.g*8e3;     % Initial condition for the geopotential field (m^2/s^2).
atm.f = 2*atm.Omega*(-atm.pts.x*sin(atm.alpha) + atm.pts.z*cos(atm.alpha)); % Coriolis force.
atm.K = 7.848e-6; % A constant (1/s)
atm.R = 4; % Wave number of the R-H wave

atm.ghm = zeros(atm.pts.nd,1);