function mt_plot(par, H)

nodes = getnodes(par);

atm.g = 9.80616;       % Gravitational constant (m/s^2).

switch par.test
  case 'galew'
    atm.gh0 = atm.g*10158.2950420088891;     % Initial condition for the geopotential field (m^2/s^2).
    zeta = mt_calc_zeta(par, H);
    ymin = 0;
  case 'tc5'
    atm.gh0 = atm.g*5960;     % Initial condition for the geopotential field (m^2/s^2).
    zeta = (H(:,4)+atm.gh0)/atm.g;
    ymin = -pi/2-0.01;
  otherwise
    error(['Unknown test case "' par.test '"']);
end

[lam,theta] = cart2sph(nodes(:,1), nodes(:,2), nodes(:,3));
tri_lt = delaunay(lam, theta);

trisurf(tri_lt,lam,theta,zeta), shading interp;
axis([-pi-0.01 pi+0.01 ymin pi/2+0.01]);
axis off;

%scatter(lam, theta, 15, zeta, 'filled')
%axis([-pi-0.01 pi+0.01 ymin pi/2+0.01]);
