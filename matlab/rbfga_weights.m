function weights = rbfga_weights(op,ep,xi,xc,opts)
%%% weights = rbfga_weights(op,ep,xi,xc,opts)
% Compute differentiation weights using the RBF-GA algorithm.
% op - operator(s) e.g., '0','x','xy' or 'L', or {'x','xy','L1'}
% ep - shape parameter
% xi - node points in 1D, 2D or 3D (n-by-d matrix, where d={1,2,3})
% xc - evaluation point (1-by-d vector)
% opts - (optional) struct containing options. [] = default
%     .verbose = [0],1 - verbose output
%     .const = [0],1   - force derivative weights to be exact for constants
%     .rtol            - tolerance for reordering [1e4*eps]
%     .m               - polynomials per order, useful when points are 
%                       non-unisolvent but reordering is unnecessary
%
% Copyright Erik Lehto, 2012-2014

[n,dim] = size(xi);
if dim>3, error('Number of dimensions must be <= 3.'); end

% Parse options
if nargin<4
    error('Too few input arguments.');
end

% Set defaults
verbose = 0;
const = 0;
rtol = 1e4*eps;
if nargin > 4
    if ~isstruct(opts), error('Variable opts must be a struct.'); end
    if isfield(opts,'verbose'), verbose = opts.verbose; end
    if isfield(opts,'const'), const = opts.const; end
    if isfield(opts,'rtol'), rtol = opts.rtol; end
    if isfield(opts,'m'), m = opts.m; end
end
if ~iscell(op), tmp = op; op = cell(1); op{1}=tmp; end % "str2cell"

% Shift and scale to unit disk
xi = xi-repmat(xc(1,:),n,1);
dzero = all(xi==0,1);
xi = xi(:,~dzero);      % Remove any zero-dimension
dim = dim-nnz(dzero);
h = max(sqrt(sum(xi.^2,2)));
xi = xi/h; ep = ep*h;

% Determine polynomial order and compute matrix
if exist('m','var')
    if sum(m) ~= n
        error('Invalid polynomial vector; sum must equal number of nodes.');
    end
    p = length(m)-1;
    P = vandermonde(xi,p);
    idx = 1:n;
else
    p = max([n-1 ceil((sqrt(n)-1.74)/0.7098) ceil((n^(1/3)-1.6)/0.555)],0);
    p = p(dim);                     % polynomial order
    P = vandermonde(xi,p);          % create multivariate polynomial matrix
    m = [ones(1,p+1); 1:p+1; cumsum(1:p+1)];
    m = m(dim,:);
    if rank(P(:,1:sum(m)))<sum(m)
        % Rank deficiency, node reordering and increase of p required
        if verbose, fprintf(1,'Reordering...   '); end
        [idx,m] = order_qr(xi,p+5,rtol);
        if verbose, fprintf(1,'Done.\n'); end
        xi = xi(idx,:);
        p = length(m)-1;
        P = vandermonde(xi,p);
    else
        idx = 1:n;
    end
    if sum(m)<n
        m = [m n-sum(m)];           % # of basis functions per order
        p = length(m)-1;
    end
end
if verbose
    disp(['# of polynomials of degree k, k=0,...: ' num2str(m,' %i')]);
end

% Total number of polynomials up to a given order
cs = cumsum([ones(1,p+1); 1:p+1; cumsum(1:p+1)],2);
cs = cs(dim,:);

% Compute constant matrices
Z = zeros(n);
for i = 1:dim
    Z = Z+xi(:,i)*xi(:,i)';
end
Z = 2*ep^2*Z;
r2 = sum(xi.^2,2);
exp_r2 = repmat(exp(-ep^2*r2'),n,1);

% Pad xi to 3D to avoid errors in computing derivatives
xi = [xi zeros(n,3-dim)];

% Create new basis, k is polynomial order of current iteration
phi = zeros(n);             % reserve memory for basis matrix
rhs = zeros(n,length(op));  % ...and for right hand side
r = 0;
for k=0:p
    rows = r+1:r+m(k+1);    % basis functions to compute this iteration
    if k==0
        B = 1;
    else
%         B = null(P(1:sum(m(1:k+1)),1:cs(k))')';
        [Q,~,~] = qr(P(1:sum(m(1:k+1)),1:cs(k))); % QR is faster than null
        B = Q(:,rows)';
    end
    nb = size(B,2);
    phi(rows,:) = exp(2*k)/ep^(2*k)*exp_r2(1:m(k+1),:).*(B*Gk(Z(1:nb,:),k));
    for i = 1:length(op)
        switch op{i}
            case '0'
                rhs(rows,i) = exp(2*k)/ep^(2*k)*...
                    (B*(Gk(0,k)*ones(nb,1)));
            case 'x'
                rhs(rows,i) = (1/h)*2*exp(2*k)/ep^(2*(k-1))*...
                    (B*(xi(1:nb,1).*Gk(0,k-1))); % d/dx
            case 'y'
                rhs(rows,i) = (1/h)*2*exp(2*k)/ep^(2*(k-1))*...
                    (B*(xi(1:nb,2).*Gk(0,k-1))); % d/dy
            case 'z'
                rhs(rows,i) = (1/h)*2*exp(2*k)/ep^(2*(k-1))*...
                    (B*(xi(1:nb,3).*Gk(0,k-1))); % d/dz
            case 'xx'
                rhs(rows,i) = (1/h^2)*2*exp(2*k)/ep^(2*(k-1))*...
                    (B*(-Gk(0,k)*ones(nb,1) + ...
                    2*ep^2*xi(1:nb,1).^2.*Gk(0,k-2))); %d2/dx2
            case 'xy'
                rhs(rows,i) = (1/h^2)*4*exp(2*k)/ep^(2*(k-2))*...
                    (B*(xi(1:nb,1).*xi(1:nb,2).*Gk(0,k-2)));
            case 'xz'
                rhs(rows,i) = (1/h^2)*4*exp(2*k)/ep^(2*(k-2))*...
                    (B*(xi(1:nb,1).*xi(1:nb,3).*Gk(0,k-2)));
            case 'yy'
                rhs(rows,i) = (1/h^2)*2*exp(2*k)/ep^(2*(k-1))*...
                    (B*(-Gk(0,k)*ones(nb,1) + ...
                    2*ep^2*xi(1:nb,2).^2.*Gk(0,k-2))); %d2/dy2
            case 'yz'
                rhs(rows,i) = (1/h^2)*4*exp(2*k)/ep^(2*(k-2))*...
                    (B*(xi(1:nb,2).*xi(1:nb,3).*Gk(0,k-2)));
            case 'zz'
                rhs(rows,i) = (1/h^2)*2*exp(2*k)/ep^(2*(k-1))*...
                    (B*(-Gk(0,k)*ones(nb,1) + ...
                    2*ep^2*xi(1:nb,3).^2.*Gk(0,k-2))); %d2/dy2
            case 'L'
                rhs(rows,i) = (1/h^2)*exp(2*k)/ep^(2*(k-1))*...
                    (B*(-2*dim*Gk(0,k)*ones(nb,1) + ...
                    4*ep^2*Gk(0,k-2)*r2(1:nb))); % Laplacian
            case {'L0','L1','L2','L3','L4','L5','L6','L7','L8','L9'}
                q = round(str2double(op{i}(2)));
                Lp = Laguerre_coeff_scaled(q, dim/2-1);
                for l = q:-1:0
                    rhs(rows,i) = rhs(rows,i) + ...
                        Lp(q-l+1)*B*((ep^2*r2(1:nb)).^l.*Gk(0,k-2*l));
                end
                rhs(rows,i) = (-4/h^2)^q*exp(2*k)*ep^(2*(q-k))*rhs(rows,i);
        end
    end
    r = r+m(k+1);
end

if verbose, disp(['Estimated condition number of interpolation matrix: ' ...
        num2str(condest(phi),' %.1e')]); end
if const
    phi = [phi ones(n,1); ones(1,n), 0];
    rhs = [rhs; zeros(1,length(op))];
end
w = phi\rhs; % Solve for differentiation weights
weights(idx,:) = w(1:n,:);

end % End of main function

function g = Gk(z,k)
if k>0
    g = exp(z).*real(gammainc(z,k));
else
    g = exp(z);
end
end

function [idx,m] = order_qr(x,p,tol)
[n,dim] = size(x);
P = vandermonde(x,p);
cs = cumsum([ones(1,p+1); 1:p+1; cumsum(1:p+1)],2);
cs = cs(dim,:);
idx = 1;
m = 1;
rem = 2:n;
i = 2;
while length(idx)<n && i<=p+1
    r = sum(m);
    % Find number of new basis functions for this stage (mj)
    [~,Rp,~] = qr(P(:,1:cs(i)));
    rankP = nnz(abs(diag(Rp))>tol);
    mj = rankP-r;
    m = [m mj];    
    for j = 1:mj
        % QR-decompose current B-matrix
        [Q,~] = qr(P(idx,1:cs(i))');
        Rnew = P(rem,1:cs(i))*Q(:,end);
        % Select new column
%         idx_new = find(abs(Rnew)>tol,1,'first');
        [~,idx_new] = max(abs(Rnew));
        idx = [idx rem(idx_new)];
        rem(idx_new) = [];
    end
    i = i+1;
end
while length(idx)<n
    [idx,m] = order_qr(x,p+3,tol);
end
end

function P = vandermonde(x,p)
[n,d] = size(x);
x = [x zeros(n,3-d)];          % pad with zeros
P = zeros(n,nchoosek(p+3,3));  % make room for 3D van der Monde
i = 1;
for k = 0:p
    for l = k:-1:0
        for m=k-l:-1:0
            P(:,i) = x(:,1).^l.*x(:,2).^m.*x(:,3).^(k-l-m);
            i = i+1;
        end
    end
end
P(:,all(P==0,1)) = [];          % truncate to correct number of columns
end

function c = Laguerre_coeff_scaled(n, alpha)
% Generalized Laguerre coefficients scaled by factorial(n)
% Code courtesy of Geert Van Damme, in a comment to File Exchange Id #15916
i = 1:n;
a = (2*i-1) + alpha;
b = sqrt(i(1:n-1).*((1:n-1) + alpha));
CM = diag(a) + diag(b,1) + diag(b,-1);
c = (-1)^n*poly(CM);
end