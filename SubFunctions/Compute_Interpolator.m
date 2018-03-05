function [N,dG,CONNECT,spElems,mspoints] = Compute_Interpolator(spCount,cellCount,x_sp,le,NN,LOC)

spElems = zeros(spCount,1);                            % index of elements where stores particles
CONNECT = zeros(spCount,4);                            % node 1=leftdown 2=righdown 3=rightup 4= leftup
N = zeros(spCount,4);                                  % Shape function
dN = zeros(spCount,8);                                 % Gradient of shape function
 
 for sp = 1:spCount
 spElems(sp) = ceil(x_sp(sp,1)/le(1))+(NN(1)-1)*(fix(x_sp(sp,2)/le(2)));   % compute vector store index elements                           
 
 CONNECT(sp,1) = spElems(sp) + floor(spElems(sp)/(NN(1)-1));
 CONNECT(sp,2) = CONNECT(sp,1)+1; 
 CONNECT(sp,3) = CONNECT(sp,2)+NN(1); 
 CONNECT(sp,4) = CONNECT(sp,1)+NN(1);

[N(sp,1),dN(sp,1),dN(sp,5)]=linearshape(x_sp(sp,1:2),LOC(CONNECT(sp,1),:),le(1,1),le(1,2));
[N(sp,2),dN(sp,2),dN(sp,6)]=linearshape(x_sp(sp,1:2),LOC(CONNECT(sp,2),:),le(1,1),le(1,2));
[N(sp,3),dN(sp,3),dN(sp,7)]=linearshape(x_sp(sp,1:2),LOC(CONNECT(sp,3),:),le(1,1),le(1,2));
[N(sp,4),dN(sp,4),dN(sp,8)]=linearshape(x_sp(sp,1:2),LOC(CONNECT(sp,4),:),le(1,1),le(1,2));

 % Build matrix of gradient of shape function of 4 nodes of the cell which
 % contain the particles
dG((sp-1)*4+[1:4],:) = [dN(sp,1) dN(sp,5);dN(sp,2) dN(sp,6);dN(sp,3) dN(sp,7);dN(sp,4) dN(sp,8)];
 end
 
 % Compute mspoints: index of particles in each element (active element)
 for c =1:cellCount
     id_sp = find(spElems==c);
     mspoints{c}=id_sp;
 end
 