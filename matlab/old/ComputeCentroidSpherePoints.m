% Author: Grady B. Wright
function xc = ComputeCentroidSpherePoints(tri,nodes)

% Local variables to make the code easier to read.
x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);

xc = 1/3*[sum(x(tri),2) sum(y(tri),2) sum(z(tri),2)];

%
% Project onto the sphere.
%

% Get the radius of the sphere the nodes lie on.
rr = sqrt(x(1)^2 + y(1)^2 + z(1)^2);

% Convert the centroids to spherical coordinates
[lam,th,r] = cart2sph(xc(:,1),xc(:,2),xc(:,3));

% Convert back to cartesian coordinates with with correct reference radius
[xc(:,1),xc(:,2),xc(:,3)] = sph2cart(lam,th,rr);
