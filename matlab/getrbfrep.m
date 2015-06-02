function [weights] = getrbfrep(par, H)

x = getnodes(par);
tree = gettree(par);
fdsize = par.fd;
ep = par.ep;

rbf = @(ep,rd2) exp(-ep^2*rd2);

N = length(x);

weights = zeros(N,fdsize);
for j=1:N
    imat = tree(:, j);

    rd2 = max(0,2*(1-x(imat,1)*x(imat,1).' - ...
                     x(imat,2)*x(imat,2).' - ...
                     x(imat,3)*x(imat,3).'));

    A(1:fdsize,1:fdsize) = rbf(ep,rd2);
    B(1:fdsize,1) = H(imat,4);
    tmpweights = A\B;

    weights(j,:) = tmpweights;
end
