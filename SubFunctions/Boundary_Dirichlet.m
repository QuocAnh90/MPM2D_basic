function [data]=Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,data)

for i=1:nfbcx
data(fbcx(i),1)     = 0;
end

for i=1:nfbcy
data(fbcy(i),2)     = 0;
end