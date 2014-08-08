function [I, J, weightsDx, weightsDy, weightsDz, weightsL] = rbfmatrix_fd_hyper_proj_ga_choose(x,ep,fdmin,fdmax,order,dim)
%%% [D,L] = rbfmatrix_fd_tree(x,ep,alpha,fdsize,order,dim)
% Requires kd-tree code.
% IN:
% x - nodes
% ep - shape parameter
% fdsize - stencil size (good choices: 31, 50, 74, 101)
% order - L = L^order
% dim - dimension of Laplacian formula
% OUT:
% DPx - sparse differentiation matrix
% DPy - sparse differentiation matrix
% DPy - sparse differentiation matrix
% L - sparse dissipation matrix

N = length(x);

weightsDx = zeros(N,fdmax);
weightsDy = zeros(N,fdmax);
weightsDz = zeros(N,fdmax);
weightsL = zeros(N,fdmax);

I = repmat((1:N)',1,fdmax);
J = knnsearch(x,x,'k',fdmax);

warning('OFF', 'MATLAB:nearlySingularMatrix')
msg_length = 0;
for j=1:N

    imat = J(j,:);
    
    xi = x(imat,:);
    P = [(x(j,2)^2+x(j,3)^2), -x(j,1)*x(j,2), -x(j,1)*x(j,3);
        -x(j,1)*x(j,2), (x(j,1)^2+x(j,3)^2), -x(j,2)*x(j,3);
        -x(j,1)*x(j,3), -x(j,2)*x(j,3), (x(j,1)^2+x(j,2)^2)];
    xp = repmat(x(j,:),fdmax,1)+xi*P;
    
    weights = zeros(fdmax,4);
    optimum = Inf*ones(1,4);
    for i=0:(fdmax-fdmin)
        n = fdmin+i;

%         in = 1:n;
        in = 2:n;
        
        opts.const = 1;
        opts.verbose = 0;

        w = rbfga_weights({'x','y','z'},ep,xp(in,:),xp(1,:),opts);
        for d=1:3
            % With center weight (in = 1:n above)
%             if sum(abs(w(1:n,d)))/sqrt(n)<optimum(d)
%                 optimum(d) = sum(abs(w(1:n,d)))/sqrt(n);
%                 weights(in,d) = w(1:n,d);

            % Without center weight (in = 2:n above)
            if sum(abs(w(1:n-1,d)))/sqrt(n)<optimum(d)
                optimum(d) = sum(abs(w(1:n-1,d)))/sqrt(n);
                weights(in,d) = w(1:n-1,d);
            end
        end
        % Center weight for L^k, differently scaled optimum
        d = 4; in = 1:n;
        w = rbfga_weights(['L' num2str(order)],ep,xp(in,:),xp(1,:),opts);
        if sum(abs(w(1:n)))<optimum(d)
            optimum(d) = sum(abs(w(1:n)));
            weights(in,d) = w(1:n);
        end
    end

    weightsDx(j,:) = weights(1:fdmax,1);
    weightsDy(j,:) = weights(1:fdmax,2);
    weightsDz(j,:) = weights(1:fdmax,3);
    weightsL(j,:) = weights(1:fdmax,4);
    
    fprintf(1,repmat('\b',1,msg_length));
    msg = [num2str(j,'%d') '/' num2str(N,'%d')];
    fprintf(1,msg);
    msg_length=numel(msg);
end
%fprintf(1,'\n');
%DPx = sparse(I,J,weightsDx,N,N); 
%DPy = sparse(I,J,weightsDy,N,N); 
%DPz = sparse(I,J,weightsDz,N,N); 
%L = sparse(I,J,weightsL,N,N);
warning('ON', 'MATLAB:nearlySingularMatrix')
