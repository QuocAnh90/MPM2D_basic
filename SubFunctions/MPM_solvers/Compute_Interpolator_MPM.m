function [N,dN,CONNECT,pElems,mpoints,NODES] = Compute_Interpolator_MPM(pCount,cellCount,x_p,le,NN,LOC)

pElems         = zeros(pCount,1);                             % index of elements where stores particles
CONNECT_TEMP    = zeros(pCount,4);                            % node 1=leftdown 2=righdown 3=rightup 4= leftup
NODES           = 4 * ones(pCount,1);                         % Number of interation nodes for each particle
N_local         = zeros(pCount,4);                            % Value of shape function
dN_local        = zeros(pCount,8);                            % Value of gradient of shape function

 for p = 1:pCount
%  pElems(p) = ceil(x_p(p,1)/le(1))+(NN(1)-1)*(fix(x_p(p,2)/le(2)));   % compute vector store index elements 
 pElems(p) = floor(x_p(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_p(p,2)/le(2)));
 
 CONNECT_TEMP(p,1) = pElems(p) + floor(pElems(p)/(NN(1)-1));
 CONNECT_TEMP(p,2) = CONNECT_TEMP(p,1)+1; 
 CONNECT_TEMP(p,3) = CONNECT_TEMP(p,2)+NN(1); 
 CONNECT_TEMP(p,4) = CONNECT_TEMP(p,1)+NN(1);
 CONNECT{p}        = [CONNECT_TEMP(p,1) CONNECT_TEMP(p,2) CONNECT_TEMP(p,3) CONNECT_TEMP(p,4)];
 
for i = 1:NODES(p)
     % Compute the shape functions and gradient of the shape functions
    [N_local(p,i),dN_local(p,i),dN_local(p,i+4)]=linearshape(x_p(p,1:2),LOC(CONNECT_TEMP(p,i),:),le(1),le(2));
    N{p}(i) = N_local(p,i);    
    dN{p}(1,i) = dN_local(p,i);
    dN{p}(2,i) = dN_local(p,i+4);
 end
 end
 
 % Compute mspoints: index of particles in each element (active element)
 for c =1:cellCount
     id_p = find(pElems==c);
     mpoints{c}=id_p;
 end
 