function mt_makecart(par)

nodes = getnodes(par.n, par.fd);
tree = gettree(par);

outname=[ '../cartfd/' num2str(par.n) '_' num2str(par.fd) '_ep' num2str(par.ep)  '.mat'];
[i,j,dx,dy,dz] = mt_rbfmatrix_fd_cart(nodes, tree, par.fd, par.ep);
save(outname, 'i','j','dx','dy','dz');
