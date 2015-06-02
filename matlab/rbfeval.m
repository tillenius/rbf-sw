function fx = rbfeval(weights,par,evalpts)

ep = par.ep;
tree = gettree(par);
x = getnodes(par);

rbf = @(ep,rd2) exp(-ep^2*rd2);

if size(evalpts,2) == 2
     [ex, ey, ez] = sph2cart(evalpts(:,1), evalpts(:,2), 1);
     evalpts = [ex, ey, ez];
end

% find closest stencil center of each evaluation point
closest = knnsearch(x,evalpts);

N = size(evalpts,1);
fx = zeros(N,1);

for j=1:N
    k = closest(j);
    imat = tree(:, k);

    % r(x) = ||x-x_k|| = sqrt((x-x_k)^2 + (y-y_k)^2 + (z-z_k)^2) = sqrt(2*(1-x'x_k))
    rd2 = max(0,2*(1-x(imat,1)*evalpts(j,1).' - ...
                     x(imat,2)*evalpts(j,2).' - ...
                     x(imat,3)*evalpts(j,3).'));

    A = rbf(ep,rd2);
    w = weights(k,:);

    fx(j) = w*A;
end
