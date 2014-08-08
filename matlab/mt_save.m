function mt_save(par)

% The design is that MATLAB reads and writes in C++ format.
% This means that indexing should be corrected to 0-based in MATLAB.

nodes = getnodes(par.n, par.fd);

infix = [num2str(par.n) ...
      '-' num2str(par.fd) ...
      '-ep' num2str(par.ep) ...
      '-o' num2str(par.order) ...
      '-gc' num2str(par.gamma_c)];

warning('off', 'MATLAB:nearlySingularMatrix')

outpath = ['../data/' par.test '-' infix];

par.dim=2;
par.n=length(nodes);
par.gamma=par.gamma_c*par.n^-par.order;

[~,~,~] = mkdir(outpath);

disp([par.test '-' infix])
switch par.test
  case 'galew'
    atm     = galew_setup(nodes);
    [uc,gh] = galew_computeInitialCondition(atm);
  case 'tc5'
    atm     = tc5_setup(nodes);
    [uc,gh] = tc5_computeInitialCondition(atm);
  otherwise
    error(['unknown test case "' par.test '"']);
end

clear nodes;

fid = fopen([outpath '/params'], 'wb');
if fid == -1
  error('Unable to open file.');
end
fwrite(fid, par.n, 'uint64');
fwrite(fid, atm.gh0, 'double');
fclose(fid);

% H contains the velocity components in x,y,z and the geopotential without
% the mean offset gh0.
H = [uc gh]';

fid = fopen([outpath '/H'], 'wb');
if fid == -1; error('Unable to open file.'); end
[m,n] = size(H);
fwrite(fid, m, 'uint64');
fwrite(fid, n, 'uint64');
fwrite(fid, H, 'double');
fclose(fid);


clear uc gh H

nd = atm.pts.nd;

tree = gettree(par);

disp('mt_rbfmatrix_fd_hyper')
[ind_i, ind_j, weightsDx, weightsDy, weightsDz, weightsL] = ...
  mt_rbfmatrix_fd_hyper([atm.pts.x atm.pts.y atm.pts.z], tree, par, atm.a);

clear tree

N = par.n;

disp('save D')

fid = fopen([outpath '/D'], 'wb');
if fid == -1; error('Unable to open file.'); end
fwrite(fid, N, 'uint64');
fwrite(fid, length(ind_i), 'uint64');
fwrite(fid, [ind_i-1 ind_j-1]', 'uint32');
fwrite(fid, [weightsDx, weightsDy, weightsDz, weightsL]', 'double');
fclose(fid);

clear weightsL


disp('gradghm')
% Gradient of the mountain.
% Compute the projected gradient of the mountain:
gradghm = zeros(nd,3);
DPx = sparse(ind_i,ind_j,weightsDx,N,N);
gradghm(:,1) = DPx*atm.ghm;
clear DPx

DPy = sparse(ind_i,ind_j,weightsDy,N,N);
gradghm(:,2) = DPy*atm.ghm;
clear DPy

DPz = sparse(ind_i,ind_j,weightsDy,N,N);
gradghm(:,3) = DPz*atm.ghm;
clear DPz

disp('save atm')
out = [atm.f atm.pts.x atm.pts.y atm.pts.z atm.pts.p_u atm.pts.p_v atm.pts.p_w atm.ghm gradghm]';

fid = fopen([outpath '/atm'], 'wb');
if fid == -1; error('Unable to open file.'); end
[m,n] = size(out);
fwrite(fid, m, 'uint64');
fwrite(fid, n, 'uint64');
fwrite(fid, out, 'double');
fclose(fid);
