function mt_plot3(par, H)

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

scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 15, H(:,4), 'filled')

switch par.test
  case 'tc5' % Draw the mountain

    % Parameters for the mountain:
    atm.lam_c = -pi/2;
    atm.thm_c = pi/6;
    atm.mR = pi/9;
    atm.hm0 = 2000;

    [atm.pts.la,atm.pts.th,r] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

    atm.pts.nd = length(H);
    atm.ghm = zeros(atm.pts.nd,1);
    r2 = (atm.pts.la-atm.lam_c).^2 + (atm.pts.th-atm.thm_c).^2;
    id = (.95*atm.mR^2 < r2) & (r2 < atm.mR^2) ;

    hold on;
    scatter3(nodes(id,1), nodes(id,2), nodes(id,3), 10, 'w', 'filled')
end

axis equal
view(50,20)
axis off
