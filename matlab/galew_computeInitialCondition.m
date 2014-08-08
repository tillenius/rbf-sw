% Initial condition for the Galewsky test case.  uc contains the velocity in Cartesian 
% coordinates with uc(:,1:3) the x,y,z direction, respectively.  gh contains
% the geopotential height field from the mean geopotential gh0.
function [uc,gh,us] = galew_computeInitialCondition(atm)

% Extract out the constants from the atm structure to make the code easier
% to read.
a = atm.a; umax = atm.umax; f = atm.f; g = atm.g; gh0 = atm.gh0;
la = atm.pts.la; th = atm.pts.th;
nd = atm.pts.nd;

ph0 = atm.ph0; ph1 = atm.ph1; ph2 = atm.ph2;
en = atm.en; alpha = atm.alpha; beta = atm.beta; hhat = atm.hhat;

us = zeros(nd,2);
ind = find(th>ph0&th<ph1); 
us(ind,1) = umax*exp(1./((th(ind)-ph0).*(th(ind)-ph1)))/en;

% Transformation for converting the latitudinal velocity to cartesian velocity.
s2c_u = [-sin(la) cos(la) zeros(size(la))];
% Transformation for converting the logitudinal velocity to cartesian velocity.
s2c_v = [-cos(la).*sin(th) -sin(la).*sin(th) cos(th)];

uc = s2c_u.*repmat(us(:,1),1,3) + s2c_v.*repmat(us(:,2),1,3);

% uu = @(th) vif(th>ph0&th<ph1,umax*exp(1./((th-ph0).*(th-ph1)))/en,0);
% rhs = @(th) a*uu(th).*(2*atm.omega*cos(th)+tan(th).*uu(th)/a);

% for i=1:nd
%     gh(i) = gh0-quadgk(rhs,-pi/2,th(i));
% end

% [ths,it] = sort(th);
% I = quadgk(rhs,-pi/2,ths(1));
% gh = zeros(nd,1);
% gh(it(1)) = gh0-I;
% for i = 2:nd
%     I = I+quadgk(rhs,ths(i-1),ths(i),'RelTol',1e-13);
%     gh(it(i)) = gh0-I;
% end
% Interpolate and evaluate the integral for the initial condition using RBFs
% gh = gh0 - rbf_integrate(rhs,[-pi/2 pi/2],th);

load int10001 % contains x and value of integral 'int' in 10001 points
gh = gh0 - interp1(x,int,th,'spline');

%ind = find(la>-pi&la<pi);
%hp(ind) = hhat*cos(th(ind)).*exp(-(la(ind)/alpha).^2).*exp(-((ph2-th(ind))/beta).^2);

hp = hhat*cos(th).*exp(-(la/alpha).^2).*exp(-((ph2-th)/beta).^2);
gh = gh + g*hp;
