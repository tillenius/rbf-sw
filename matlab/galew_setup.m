%
% Set up for the Galewsky test case.
%
function atm = galew_setup(nfile)

if ischar(nfile)
   nodes = load(nfile);
else
   nodes = nfile;
end
[atm.pts.la,atm.pts.th] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
[nodes(:,1),nodes(:,2),nodes(:,3)] = sph2cart(atm.pts.la,atm.pts.th,ones(size(atm.pts.la)));

atm.pts.nodes = nodes(:,1:3);
atm.pts.x = nodes(:,1); atm.pts.y = nodes(:,2); atm.pts.z = nodes(:,3);
atm.pts.nd = size(atm.pts.x,1);
if size(nodes,2) == 4   % Extract the quadrature weights.
   atm.pts.wghts = nodes(:,4)';
else
   atm.pts.wghts = repmat(atm.pts.nd/(4*pi),[1 atm.pts.nd]);
end

% Variables for projecting an arbitrary Cartesian vector onto the surface
% of the sphere.
x2 = nodes(:,1).^2; xy = nodes(:,1).*nodes(:,2);
y2 = nodes(:,2).^2; xz = nodes(:,1).*nodes(:,3);
z2 = nodes(:,3).^2; yz = nodes(:,2).*nodes(:,3);
atm.pts.p_u = [1-x2  -xy   -xz];
atm.pts.p_v = [-xy  1-y2   -yz];
atm.pts.p_w = [-xz   -yz  1-z2];

atm.umax = 80; % Speed of rotation in meters/second
atm.a = 6.37122e6;     % Mean radius of the earth (meters).
atm.omega = 7.292e-5;  % Rotation rate of the earth (1/seconds).
atm.g = 9.80616;       % Gravitational constant (m/s^2).
atm.gh0 = atm.g*10158.2950420088891;     % Initial condition for the geopotential field (m^2/s^2).
atm.f = 2*atm.omega*atm.pts.z; % Coriolis force.

atm.ph0 = pi/7;
atm.ph1 = pi/2 - atm.ph0;
atm.en = exp(-4/(atm.ph1-atm.ph0)^2);

atm.alpha = 1/3;
atm.beta = 1/15;
atm.ph2 = pi/4;
atm.hhat = 120;

atm.ghm = zeros(atm.pts.nd,1);