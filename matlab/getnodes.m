function nodes=getnodes(nodes, fdsize)

if isstruct(nodes)
    fdsize = nodes.fd;
    nodes = nodes.n;
end

if isscalar(nodes)
    nodes = ['../nodesets/' num2str(nodes) '_' num2str(fdsize) '.mat'];
end

if ischar(nodes)
    nodes=load(nodes);
    if isstruct(nodes)
        nodes=nodes.nodes;
    end
end
