% %% load and compute error vs. DWD, local or global RBFs
% nvals = (60:10:160).^2;
nvals = [3600        6400       10000       16900       25600];

err_fd17 = zeros(length(nvals),3);
err_e_fd17 = zeros(length(nvals),1);
err_m_fd17 = zeros(length(nvals),1);

err_fd31 = zeros(length(nvals),3);
err_e_fd31 = zeros(length(nvals),1);
err_m_fd31 = zeros(length(nvals),1);

err_fd50 = zeros(length(nvals),3);
err_e_fd50 = zeros(length(nvals),1);
err_m_fd50 = zeros(length(nvals),1);

err_fd101 = zeros(length(nvals),3);
err_e_fd101 = zeros(length(nvals),1);
err_m_fd101 = zeros(length(nvals),1);

refname = 'data/TC5_n163842_ep14.0671_fd31';
% refname = 'data/TC5_n163842_ep14.0671_fd31g';
% refname = 'data\TC5_n25600_ep9.12_fd101g';
load([refname '.mat']);
xe = atm.pts.x; ye = atm.pts.y; ze = atm.pts.z;
h = (H(:,4)+atm.gh0)/atm.g;
% wghts = atm.pts.wghts;
trie = TriangulateSpherePoints(atm.pts.nodes);
wghts = ComputeTriQuadWeights(trie,atm.pts.nodes);
trisa = ComputeTriSurfaceArea(trie,atm.pts.nodes);

for i=1:length(nvals)
    n = nvals(i);

    fd = 17;
    ep = 0.026*sqrt(n)-0.08;
    bsname = ['../tc5/data/TC5_n' num2str(n) '_ep' num2str(ep) ...
        '_fd' num2str(fd)];    load([bsname '.mat']);
    x = atm.pts.x; y = atm.pts.y; z = atm.pts.z;
    he = rbf_interp_fd([x y z],ep,fd,(H(:,4)+atm.gh0)/atm.g,[xe ye ze]);
    diff = he - h;
    errs = [sum(trisa.*sum(wghts.*abs(diff(trie)),2)) ...
        sqrt(sum(trisa.*sum(wghts.*diff(trie).^2,2))) norm(diff,inf)]./...
        [sum(trisa.*sum(wghts.*abs(h(trie)),2)) ...
        sqrt(sum(trisa.*sum(wghts.*h(trie).^2,2))) norm(h,inf)];
    err_fd17(i,:) = errs;
    err_e_fd17(i) = max(abs(erre));
    err_m_fd17(i) = max(abs(errm));
    disp(num2str([n fd err_fd17(i,2) err_m_fd17(i) err_e_fd17(i)]));

    
    ep=0.035*sqrt(n)-0.1;
    fd = 31;
%     bsname = ['data/TC5_n' num2str(n) '_ep' num2str(ep) ...
%         '_fd' num2str(fd) 'g'];
    bsname = ['data/TC5_n' num2str(n) '_ep' num2str(ep) ...
        '_fd' num2str(fd)];    load([bsname '.mat']);
    x = atm.pts.x; y = atm.pts.y; z = atm.pts.z;
    he = rbf_interp_fd([x y z],ep,fd,(H(:,4)+atm.gh0)/atm.g,[xe ye ze]);
    diff = he - h;
    errs = [sum(trisa.*sum(wghts.*abs(diff(trie)),2)) ...
        sqrt(sum(trisa.*sum(wghts.*diff(trie).^2,2))) norm(diff,inf)]./...
        [sum(trisa.*sum(wghts.*abs(h(trie)),2)) ...
        sqrt(sum(trisa.*sum(wghts.*h(trie).^2,2))) norm(h,inf)];
%     errs = [wghts*abs(diff) sqrt(wghts*diff.^2) norm(diff,inf)]./...
%         [wghts*abs(h) sqrt(wghts*h.^2) norm(h,inf)];
    err_fd31(i,:) = errs;
    err_e_fd31(i) = max(abs(erre));
    err_m_fd31(i) = max(abs(errm));
    disp(num2str([n fd err_fd31(i,2) err_m_fd31(i) err_e_fd31(i)]));
    
    fd = 50;
    ep = 0.044*sqrt(n)-0.14;
    bsname = ['../tc5/data/TC5_n' num2str(n) '_ep' num2str(ep) ...
        '_fd' num2str(fd)];    load([bsname '.mat']);
    x = atm.pts.x; y = atm.pts.y; z = atm.pts.z;
    he = rbf_interp_fd([x y z],ep,fd,(H(:,4)+atm.gh0)/atm.g,[xe ye ze]);
    diff = he - h;
    errs = [sum(trisa.*sum(wghts.*abs(diff(trie)),2)) ...
        sqrt(sum(trisa.*sum(wghts.*diff(trie).^2,2))) norm(diff,inf)]./...
        [sum(trisa.*sum(wghts.*abs(h(trie)),2)) ...
        sqrt(sum(trisa.*sum(wghts.*h(trie).^2,2))) norm(h,inf)];
    err_fd50(i,:) = errs;
    err_e_fd50(i) = max(abs(erre));
    err_m_fd50(i) = max(abs(errm));
    disp(num2str([n fd err_fd50(i,2) err_m_fd50(i) err_e_fd50(i)]));
    
    ep = 0.058*sqrt(n)-0.16;
    
    fd = 101;
%     bsname = ['data/TC5_n' num2str(n) '_ep' num2str(ep) ...
%         '_fd' num2str(fd) 'g'];
    bsname = ['data/TC5_n' num2str(n) '_ep' num2str(ep) ...
        '_fd' num2str(fd)];
    load([bsname '.mat']);
    x = atm.pts.x; y = atm.pts.y; z = atm.pts.z;
    he = rbf_interp_fd([x y z],ep,fd,(H(:,4)+atm.gh0)/atm.g,[xe ye ze]);
    diff = he - h;
    errs = [sum(trisa.*sum(wghts.*abs(diff(trie)),2)) ...
        sqrt(sum(trisa.*sum(wghts.*diff(trie).^2,2))) norm(diff,inf)]./...
        [sum(trisa.*sum(wghts.*abs(h(trie)),2)) ...
        sqrt(sum(trisa.*sum(wghts.*h(trie).^2,2))) norm(h,inf)];
%     errs = [wghts*abs(diff) sqrt(wghts*diff.^2) norm(diff,inf)]./...
%         [wghts*abs(h) sqrt(wghts*h.^2) norm(h,inf)];
    err_fd101(i,:) = errs;
    err_e_fd101(i) = max(abs(erre));
    err_m_fd101(i) = max(abs(errm));
    
    disp(num2str([n fd err_fd101(i,2) err_m_fd101(i) err_e_fd101(i)]));
       
%     disp(['fd, i=' num2str(i) ',l2=' num2str(errs(2)) '.']);

end