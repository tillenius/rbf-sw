% Initial condition for test case 6.  uc contains the velocity in Cartesian 
% coordinates with uc(:,1:3) the x,y,z direction, respectively.  gh contains
% the geopotential height field from the mean geopotential gh0.
function [uc,gh,us] = tc6_computeInitialCondition(atm)

% Extract out the constants from the atm structure to make the code easier
% to read.
a = atm.a; omega = atm.omega; Omega = atm.Omega;
K = atm.K; R = atm.R; gh0 = atm.gh0;
la = atm.pts.la; th = atm.pts.th;

us = [a*omega*cos(th)+a*K*cos(th).^(R-1).*(R*sin(th).^2-cos(th).^2).*cos(R*la) ...
    -a*K*R*cos(th).^(R-1).*sin(th).*sin(R*la)];

% Transformation for converting the latitudinal velocity to cartesian velocity.
s2c_u = [-sin(la) cos(la) zeros(size(la))];
% Transformation for converting the logitudinal velocity to cartesian velocity.
s2c_v = [-cos(la).*sin(th) -sin(la).*sin(th) cos(th)];
uc = s2c_u.*repmat(us(:,1),1,3) + s2c_v.*repmat(us(:,2),1,3);

A = 0.5*omega*(2*Omega+omega)*cos(th).^2 + 0.25*K^2*cos(th).^(2*R).*...
    ((R+1)*cos(th).^2 + (2*R^2-R-2) - 2*R^2*cos(th).^(-2));
B = (2*(Omega+omega)*K)*cos(th).^R.*((R^2+2*R+2)-(R+1)^2*cos(th).^2)/((R+1)*(R+2));
C = 0.25*K^2*cos(th).^(2*R).*((R+1)*cos(th).^2 - (R+2));

gh = a^2*A + a^2*B.*cos(R*la) + a^2*C.*cos(2*R*la); % height field without mean offset