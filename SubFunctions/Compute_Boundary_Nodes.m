function [nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC)
% Input
% nodeCount: total number of nodes
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction

% Output
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction


nfbcx                   = 0              ;  % initial number of fixed nodes in x direction
nfbcy                   = 0              ;  % initial number of fixed nodes in y direction
fbcx = []; fbcy = [];                       % vector store the index of boundary nodes

for n=1:nodeCount
    if LOC(n,1)<=0.05 || LOC(n,1)>=11
        nfbcx = nfbcx+1;
        fbcx = [fbcx n];
    end
end

for n=1:nodeCount
    if LOC(n,2)<=0 || LOC(n,2)>=0.15
        nfbcy = nfbcy+1;
        fbcy = [fbcy n];
    end
end