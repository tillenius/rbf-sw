% Author: Grady B. Wright
function [wghts,gwghts] = ComputeTriQuadWeights(tri,nodes)

xc = ComputeCentroidSpherePoints(tri,nodes);

% Local variables to make the code easier to read.
x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);

% For each triangle compute the distance the centroid is each vertex.
for j=1:3
   dp = (x(tri(:,j)).*xc(:,1) + y(tri(:,j)).*xc(:,2) + z(tri(:,j)).*xc(:,3));
   % Weights based on Euclidean distance.
   wghts(:,j) = sqrt( 2*(1-dp) );
   % Weights based on Geodesic distance.
   gwghts(:,j) = acos(dp);
end
wghts = wghts./repmat(sum(wghts,2),[1 3]);
gwghts = gwghts./repmat(sum(gwghts,2),[1 3]);


