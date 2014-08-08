function tree=gettree(par)

treename = ['../tree/' num2str(par.n) '_' num2str(par.fd) '.mat'];
tree = load(treename);
if isstruct(tree)
    tree = tree.neighbors;
end

