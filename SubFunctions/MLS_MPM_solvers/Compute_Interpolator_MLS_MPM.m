function [N,dN,N_MLS,dN_MLS,N_g,N_c,dN_c,CONNECT,CONNECT_G,CONNECT_C,spElems,mspoints,NODES,GAUSS,NODES_G] = Compute_Interpolator_MLS_MPM(spCount,cellCount,nodeCount,x_sp,le,NN,LOC,LOCC)

%% Output
% Mapping from particles to nodes using linear basis N{p}(i); dN{p}(1:2,i)
% Mapping from particles to nodes using moving least square N_MLS{p}(i)
% Mapping from particles to centroids using quadratic Bspline Ng{p}(c)
% Mapping from centroids to nodes using linear basis Nc{p}(i); dNc{p}(1:2,i)


%% Parameter generations
spElems         = zeros(spCount,1);              % index of elements where stores particles
CONNECT_TEMP    = zeros(spCount,4);              % node 1=leftdown 2=righdown 3=rightup 4= leftup
CONNECT_TEMP_G  = zeros(spCount,9);              % gauss number from left to right from bottom to top
CONNECT_TEMP_C  = zeros(cellCount,4);
NODES           = 4 * ones(spCount,1);           % Number of interation nodes for each particle
GAUSS           = 9 * ones(spCount,1);
cElems          = zeros(cellCount,1);
NODES_G         = 4 * ones(cellCount,1);          

% Funtion for mapping from particles to nodes using linear basis
N_local         = zeros(spCount,4);              % Value of shape function
dNx_local       = zeros(spCount,4);              % Value of gradient of shape function
dNy_local       = zeros(spCount,4); 

% Funtion for mapping from particles to nodes using moving least square
% Least square matrix
A               = cell(nodeCount,1);
dAx             = cell(nodeCount,1);
dAy             = cell(nodeCount,1);

% Function for mapping from particles to centroids using quadratic Bspline
Ng_local         = zeros(spCount,9);              % Value of shape function
% Least square matrix
Ac               = cell(cellCount,1);

% Funtion for mapping from centroids to nodes using linear basis
Nc_local         = zeros(cellCount,4);              % Value of shape function
dNcx_local       = zeros(cellCount,4);              % Value of gradient of shape function
dNcy_local       = zeros(cellCount,4); 

% Order of MLS
r = 3;
for i = 1:nodeCount
    A{i} = zeros(r,r);
    dAx{i} = zeros(r,r);
    dAy{i} = zeros(r,r);
end
 
%% Compute the index of the node for the particles 
% and the index of the centroids to particles

 for p = 1:spCount
%  spElems(p) = ceil(x_sp(p,1)/le(1))+(NN(1)-1)*(fix(x_sp(p,2)/le(2)));   % compute vector store index elements 
 
 % Node index for particle
 CONNECT_TEMP(p,1) = spElems(p) + ceil(spElems(p)/(NN(1)-1)) - 1;
 CONNECT_TEMP(p,2) = CONNECT_TEMP(p,1)+1; 
 CONNECT_TEMP(p,3) = CONNECT_TEMP(p,2)+NN(1); 
 CONNECT_TEMP(p,4) = CONNECT_TEMP(p,1)+NN(1);
 CONNECT{p}        = [CONNECT_TEMP(p,1) CONNECT_TEMP(p,2) CONNECT_TEMP(p,3) CONNECT_TEMP(p,4)];
 
 % Centroid index for particle
 CONNECT_TEMP_G(p,1) =  spElems(p)-(NN(1)-1)-1;
 CONNECT_TEMP_G(p,2) =  spElems(p)-(NN(1)-1);
 CONNECT_TEMP_G(p,3) =  spElems(p)-(NN(1)-1)+1;
 CONNECT_TEMP_G(p,4) =  spElems(p)-1;
 CONNECT_TEMP_G(p,5) =  spElems(p);
 CONNECT_TEMP_G(p,6) =  spElems(p)+1;
 CONNECT_TEMP_G(p,7) =  spElems(p)+(NN(1)-1)-1;
 CONNECT_TEMP_G(p,8) =  spElems(p)+(NN(1)-1);
 CONNECT_TEMP_G(p,9) =  spElems(p)+(NN(1)-1)+1;
 CONNECT_G{p} = [CONNECT_TEMP_G(p,1) CONNECT_TEMP_G(p,2) CONNECT_TEMP_G(p,3) CONNECT_TEMP_G(p,4) CONNECT_TEMP_G(p,5) CONNECT_TEMP_G(p,6) CONNECT_TEMP_G(p,7) CONNECT_TEMP_G(p,8) CONNECT_TEMP_G(p,9)];
 end
 
 % Node index for centroid
 for c = 1:cellCount
 cElems(c) = ceil(LOCC(c,1)/le(1))+(NN(1)-1)*(fix(LOCC(c,2)/le(2)));
 CONNECT_TEMP_C(c,1) = cElems(c) + ceil(cElems(c)/(NN(1)-1)) - 1;
 CONNECT_TEMP_C(c,2) = CONNECT_TEMP_C(c,1)+1; 
 CONNECT_TEMP_C(c,3) = CONNECT_TEMP_C(c,2)+NN(1); 
 CONNECT_TEMP_C(c,4) = CONNECT_TEMP_C(c,1)+NN(1);
 CONNECT_C{c}        = [CONNECT_TEMP_C(c,1) CONNECT_TEMP_C(c,2) CONNECT_TEMP_C(c,3) CONNECT_TEMP_C(c,4)];
 end
     
%% Compute mspoints: index of particles in each element (active element)
 for c =1:cellCount
     id_p = find(spElems==c);
     mspoints{c}=id_p;
 end
 
%% Linear basis function from particle to nodes
 for p = 1:spCount
for i = 1:NODES(p)
    npid = CONNECT_TEMP(p,i);
     % Compute the shape functions and gradient of the shape functions
    [N_local(p,i),dNx_local(p,i),dNy_local(p,i)]=linearshape(x_sp(p,1:2),LOC(npid,:),le(1),le(2));
    N{p}(i) = N_local(p,i);    
    dN{p}(1,i) = dNx_local(p,i);
    dN{p}(2,i) = dNy_local(p,i);
end
end
 
%% Least square function from particle to nodes
% Compute A
for p = 1:spCount
     if r==1
         Pxy = [1];
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];    
     end
         
     PPxy = Pxy * Pxy';
     
     for i = 1:NODES(p)
         npid = CONNECT{p}(i);
         A{npid} = A{npid} + N_local(p,i) * PPxy;
         dAx{npid} = dAx{npid} + dNx_local(p,i) * PPxy;
         dAy{npid} = dAy{npid} + dNy_local(p,i) * PPxy;    
     end
end
 
% Compute shape function and gradients
 for p = 1:spCount
     if r==1
         Pxy = [1];
         dpdx = [0];
         dpdy = [0];
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];   
         dpdx = [0 ; 1 ; 0];
         dpdy = [0 ; 0 ; 1];
     end
     
     for i = 1:NODES(p)
         npid = CONNECT{p}(i);
         if r==1
         pn = [1];
         elseif r==3
         pn = [1 ; LOC(npid,1) ; LOC(npid,2)];   
         end
         
         rn = A{npid} \ pn;
         N_MLS{p}(i) = rn' * N_local(p,i) * Pxy;
         
         drnx = A{npid} \ (dpdx - dAx{npid} * rn);
         dN_MLS{p}(1,i) = drnx' * N_local(p,i) * Pxy +  rn' * dNx_local(p,i) * Pxy;
         
         drny = A{npid} \ (dpdy - dAy{npid} * rn);
         dN_MLS{p}(2,i) = drny' * N_local(p,i) * Pxy +  rn' * dNy_local(p,i) * Pxy;        
     end
 end
 
%% Quadratic-Bspline for mapping from particle to centroids
 for p = 1:spCount
for c = 1:GAUSS(p)
    gpid = CONNECT_TEMP_G(p,c);
     % Compute the shape functions and gradient of the shape functions
%     [Ng_local(p,c)]=linearshape(x_sp(p,1:2),LOCC(gpid,:),le(1),le(2));
    [Ng_local(p,c)]=Quadratic_Bspline(x_sp(p,1:2),LOCC(gpid,:),le(1),le(2));
    N_g{p}(c) = Ng_local(p,c);    
end
 end
 
% Least square function from particle to centroids
% Compute A
for p = 1:spCount
     if r==1
         Pxy = 1;
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];    
     elseif r==6
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
     end
         
     PPxy = Pxy * Pxy';
     
     for c = 1:GAUSS(p)
         cpid = CONNECT_TEMP_G(p,c);
         Ac{cpid} = Ac{cpid} + Ng_local(p,c) * PPxy;
     end
end

% Compute shape function and gradients
 for p = 1:spCount
     if r==1
         Pxy = 1;
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];   
     elseif r==6
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
     end
     
     for c = 1:GAUSS(p)
         cpid = CONNECT_TEMP_G(p,c);
                
         if r==1
         pc = [1];
         elseif r==3
         pc = [1 ; LOCC(cpid,1) ; LOCC(cpid,2)];   
         elseif r==6
         pc = [1 ; LOCC(cpid,1) ; LOCC(cpid,2) ; LOCC(cpid,1)^2 ; LOCC(cpid,1)*LOCC(cpid,2) ; LOCC(cpid,2)^2];   
         end
         
         if r_c(cpid)==6
         
         rc = Ac{cpid} \ pc;
         N_g{p}(c) = rc' * Ng_local(p,c) * Pxy;  
         
         elseif r_c(cpid)==3
         rc = Ac{cpid}(1:3,1:3) \ pc(1:3);
         N_g{p}(c) = rc' * Ng_local(p,c) * Pxy(1:3);  
         end
         
     end
 end
 
 %% Linear for mapping from centroids to node
 for c = 1:cellCount
     for i=1:NODES_G(c)
     npid = CONNECT_TEMP_C(c,i);
     % Compute the shape functions and gradient of the shape functions
    [Nc_local(c,i),dNcx_local(c,i),dNcy_local(c,i)]=linearshape(LOCC(c,:),LOC(npid,:),le(1),le(2));
    N_c{c}(i) = Nc_local(c,i);    
    dN_c{c}(1,i) = dNcx_local(c,i);
    dN_c{c}(2,i) = dNcy_local(c,i);
     end
 end
 
 