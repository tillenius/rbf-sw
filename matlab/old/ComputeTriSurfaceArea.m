% Author: Grady B. Wright
function sa = ComputeTriSurfaceArea(tri,nodes)

% Local variables to make the code easier to read.
x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);

% 
% Local variables for the coordinate of each triangle.
%
xt = x(tri);
yt = y(tri);
zt = z(tri);

clear x y z;

%
% Calculate the square of the radius for the sphere defined by each triangle.
%
R2 = xt(:,1).^2 + yt(:,1).^2 + zt(:,1).^2;
R = sqrt(R2);

% 
% Now for each sphere defined by the triangle compute its area.  Note that
% it may be usefull to look at the web page http://mathworld.wolfram.com/SphericalTrigonometry.html
% to understand all the mathematics and the definitions for what follows in
% the calculations below.
%
c = real(acos((xt(:,1).*xt(:,2) + yt(:,1).*yt(:,2) + zt(:,1).*zt(:,2))./R2)).*R;
b = real(acos((xt(:,1).*xt(:,3) + yt(:,1).*yt(:,3) + zt(:,1).*zt(:,3))./R2)).*R;
a = real(acos((xt(:,2).*xt(:,3) + yt(:,2).*yt(:,3) + zt(:,2).*zt(:,3))./R2)).*R;

% Calculate the angles of each of the spherical triangles.
A = acos((cos(a) - cos(b).*cos(c))./(sin(b).*sin(c)));
B = acos((cos(b) - cos(c).*cos(a))./(sin(c).*sin(a)));
C = acos((cos(c) - cos(a).*cos(b))./(sin(a).*sin(b)));

% Calculate the area of the triangles (see http://mathworld.wolfram.com/SphericalTriangle.html)
sa = R2.*(A + B + C - pi);


