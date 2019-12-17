function [N,dN,CONNECT,pElems,mpoints,NODES] = Compute_Interpolator_BSMPM(pCount,cellCount,x_sp,le,NN,LOC)

pElems          = zeros(pCount,1);                            % index of elements where stores particles
CONNECT_TEMPX   = zeros(pCount,3);
CONNECT_TEMPY   = zeros(pCount,3); 
NODES           = 3 * ones(pCount,1);                         % Number of interation nodes for each particle

% N_localx         = zeros(pCount,3);                            % Value of shape function
% dN_localx        = zeros(pCount,3);                            % Value of gradient of shape function
% N_localy         = zeros(pCount,3);                            % Value of shape function
% dN_localy        = zeros(pCount,3);                            % Value of gradient of shape function


% Knot vector reconstruction
     deg = 2;
     
     XiCount = 1/le(1)+5; 
     Xi = zeros(XiCount,1);
     Xi(1:3,1) = 2*le(1);
     for i=1:1/le(1)
     Xi(3+i,1) =  Xi(2+i,1)+le(1);
     end
     Xi(XiCount-1:XiCount) = Xi(XiCount-2);
     
     YiCount = 3;
     Yi = zeros(YiCount,1);
     Yi(1:3,1) = 1*le(2);
     Yi(4:6,1) = 2*le(2);
     
 for p = 1:pCount
%  pElems(p) = ceil(x_p(p,1)/le(1))+(NN(1)-1)*(fix(x_p(p,2)/le(2)));   % compute vector store index elements 
 pElems(p) = floor(x_sp(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_sp(p,2)/le(2)));
 
 % Knot vector index
 CONNECT_TEMPX(p,1) = floor(x_sp(p,1)/le(1)+1) - 2;
 CONNECT_TEMPX(p,2) = floor(x_sp(p,1)/le(1)+1) - 1;
 CONNECT_TEMPX(p,3) = floor(x_sp(p,1)/le(1)+1) - 0;
 CONNECTX{p}        = [CONNECT_TEMPX(p,1) CONNECT_TEMPX(p,2) CONNECT_TEMPX(p,3)];
 
 CONNECT_TEMPY(p,1) = floor(x_sp(p,2)/le(2)+1) - 1;
 CONNECT_TEMPY(p,2) = floor(x_sp(p,2)/le(2)+1) - 0;
 CONNECT_TEMPY(p,3) = floor(x_sp(p,2)/le(2)+1) + 1;
 CONNECTY{p}        = [CONNECT_TEMPY(p,1) CONNECT_TEMPY(p,2) CONNECT_TEMPY(p,3)];
 
     [N_localx , dN_localx] = test_new_bspline(Xi,deg,x_sp(p,1));
     
     for i = 1:NODES(p)
         Nx{p}(i) = N_localx(CONNECT_TEMPX(p,i));
         dNx{p}(1,i) = dN_localx(CONNECT_TEMPX(p,i));
     end
     
     [N_localy , dN_localy] = test_new_bspline(Yi,deg,x_sp(p,2));
     
     for i = 1:NODES(p)
         Ny{p}(i) = N_localy(CONNECT_TEMPY(p,i));
         dNy{p}(1,i) = dN_localy(CONNECT_TEMPY(p,i));
     end
 end
 
 
%  Compute mspoints: index of particles in each element (active element)
 for c =1:cellCount
     id_p = find(pElems==c);
     mpoints{c}=id_p;
 end