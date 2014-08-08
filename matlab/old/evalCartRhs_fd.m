% Evaluates the RHS (spatial derivatives) for the Cartesian RBF formulation of 
% the shallow water equations with projected gradients.
% This function applies to Test Case 5, which contains bottom topography
function F = evalCartRhs_fd(~,H,DPx,DPy,DPz,L,atm,gradghm, ind,weights)

nd = atm.pts.nd;
H = reshape(H,nd,4);
% Extract out some constants from the atm structure to make the code
% easier to read.
x = atm.pts.x; y = atm.pts.y; z = atm.pts.z;
f = atm.f; a = atm.a;  gh0 = atm.gh0;

% Compute the (projected) Cartesian derivatives applied to the velocity
% and geopotential.
Tx = DPx*H/a;
Ty = DPy*H/a;
Tz = DPz*H/a;

%
% This is the computation for the right hand side of the (Cartesian) 
% momentum equations.
%
p = -(H(:,1).*Tx(:,1) + H(:,2).*Ty(:,1) + H(:,3).*Tz(:,1) + (f).*(y.*H(:,3) - z.*H(:,2)) + Tx(:,4));
q = -(H(:,1).*Tx(:,2) + H(:,2).*Ty(:,2) + H(:,3).*Tz(:,2) + (f).*(z.*H(:,1) - x.*H(:,3)) + Ty(:,4));
s = -(H(:,1).*Tx(:,3) + H(:,2).*Ty(:,3) + H(:,3).*Tz(:,3) + (f).*(x.*H(:,2) - y.*H(:,1)) + Tz(:,4));

% Project the momentum equations onto the surface of the sphere.
F(:,1) = atm.pts.p_u(:,1).*p + atm.pts.p_u(:,2).*q + atm.pts.p_u(:,3).*s;
F(:,2) = atm.pts.p_v(:,1).*p + atm.pts.p_v(:,2).*q + atm.pts.p_v(:,3).*s;
F(:,3) = atm.pts.p_w(:,1).*p + atm.pts.p_w(:,2).*q + atm.pts.p_w(:,3).*s;

% Right-hand side for the geopotential (Does not need to be projected, this
% has already been accounted for in the DPx, DPy, and DPz operators for
% this equation).
F(:,4) = -(  H(:,1).*(Tx(:,4) - gradghm(:,1)) ...
           + H(:,2).*(Ty(:,4) - gradghm(:,2)) ...
           + H(:,3).*(Tz(:,4) - gradghm(:,3)) ...
           + (H(:,4)+gh0-atm.ghm).*(Tx(:,1) + Ty(:,2) + Tz(:,3)));

% Apply the hyper-viscosity, either once or twice.
F = F + L*H;
% F = F - L*(L*H);
F = reshape(F,4*nd,1);
