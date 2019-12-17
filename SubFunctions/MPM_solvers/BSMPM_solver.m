function[v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = BSMPM_solver(CModel,CModel_parameter,...
    nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
    nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt)

%% Store particles into cell
% N: Shape function
% dN: Gradient of shape function
% spElems(p): element index where "p" locates
% mspoints(e): all particles indexes where locate in the cell "e"

spElems         = zeros(spCount,1);                            % index of elements where stores particles
CONNECT_TEMPX   = zeros(spCount,3);
CONNECT_TEMPY   = zeros(spCount,3); 

% Knot vector reconstruction
     deg = 2;
     Nelemsx = (NN(1)) - 4 + 1;
     Nelemsy = 3;
     
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
     
     for p = 1:spCount
     % compute vector store index elements 
     spElems(p) = floor(x_sp(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_sp(p,2)/le(2)));

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
     
     for i = 1:3
         Nx{p}(i) = N_localx(CONNECT_TEMPX(p,i));
         dNx{p}(1,i) = dN_localx(CONNECT_TEMPX(p,i));
     end
     
     [N_localy , dN_localy] = test_new_bspline(Yi,deg,x_sp(p,2));
     
     for i = 1:3
         Ny{p}(i) = N_localy(CONNECT_TEMPY(p,i));
         dNy{p}(1,i) = dN_localy(CONNECT_TEMPY(p,i));
     end
     end
 
%% Mapping from particle to nodes
% nmass_si: nodal mass
% nmomentum_si: nodal momentum
% niforce_si: nodal internal force
% neforce_si: nodal external force
% traction_si: nodal traction

% Interpolation from particle to grid task
% Node variables
nmass                = cell(Nelemsx,Nelemsy);                   % Nodal Mass
nmomentum            = cell(Nelemsx,Nelemsy);                   % Nodal Momentum
niforce              = cell(Nelemsx,Nelemsy);                   % Nodal Internal force
neforce              = cell(Nelemsx,Nelemsy);                   % Nodal External force
traction             = cell(Nelemsx,Nelemsy);                   % Nodal Traction
nforce               = cell(Nelemsx,Nelemsy);                   % Total force
nvelo                = cell(Nelemsx,Nelemsy);                   % Nodal Velocity

% Nelemsx: number of basis function in x direction
% Nelemsy: number of basis function in y direction

% Generate the nodal data
for i = 1:Nelemsx
    for j = 1:Nelemsy
        nmass{i,j}                = 0;                   
        nmomentum{i,j}            = zeros(1,2);
        niforce{i,j}              = zeros(1,2);  
        neforce{i,j}              = zeros(1,2);    
        traction{i,j}             = zeros(1,2);      
    end
end

% Interpolation from particle to node
 for p=1:spCount
 % Build stress tensor
 SSP = [s_sp(p,1) s_sp(p,3);s_sp(p,3) s_sp(p,2)];
 
 for i = 1:3
     xpid = CONNECTX{p}(i); % global x coordinate of basis function
     for j = 1:3
         ypid = CONNECTY{p}(j); % global y coordinate of basis function
         
         % Mass
         nmass{xpid,ypid}       = nmass{xpid,ypid} + m_sp(p) * Nx{p}(i) * Ny{p}(j);
         
         % Momentum
         nmomentum{xpid,ypid}   = nmomentum{xpid,ypid} + m_sp(p) * v_ssp(p,:) * Nx{p}(i) * Ny{p}(j);
         
         % Internal force
         niforce{xpid,ypid}     = niforce{xpid,ypid} - V_sp(p) * (SSP * [dNx{p}(i)*Ny{p}(j);Nx{p}(i)*dNy{p}(j)])';
         
         % External force
         neforce{xpid,ypid}     = neforce{xpid,ypid} + b_sp(p,:) * m_sp(p) * Nx{p}(i) * Ny{p}(j); 
         
         % Traction
         traction{xpid,ypid}(1) = traction{xpid,ypid}(1) + V_sp(p) * ptraction_sp(p,1) * Nx{p}(i) * Ny{p}(j) / le(1);
         traction{xpid,ypid}(2) = traction{xpid,ypid}(2) + V_sp(p) * ptraction_sp(p,2) * Nx{p}(i) * Ny{p}(j) / le(2);
     end
 end
 end

 % Check conservation of mass
 Mass = 0;
 Momentum = 0;
 for i=1:Nelemsx
     for j = 1:Nelemsy
         Mass = Mass + nmass{i,j};
         Momentum = Momentum + nmomentum{i,j};
         nvelo{i,j} = nmomentum{i,j} / nmass{i,j}; 
     end
 end

  test = zeros(Nelemsx,1);
 for i = 1:Nelemsx
     test(i) = nvelo{i,2}(1);
 end
 
 
%% Update momentum
% Update force and momentum
for i = 1:Nelemsx
    for j = 1:Nelemsy
        nforce{i,j}     	= niforce{i,j} + neforce{i,j} + traction{i,j};
    end
end

for i = 1:Nelemsx
    for j = 1:Nelemsy
        nmomentum{i,j}      = nmomentum{i,j} + nforce{i,j}*dt;
        
        % Boundary condition
        if i == 1 || i == Nelemsx
            nforce{i,j}(1)      = 0;
            nmomentum{i,j}(1)   = 0;
        end
        
        if j == 1 || j == Nelemsy
            nforce{i,j}(2)      = 0;
            nmomentum{i,j}(2)   = 0;
        end
    end
end

%% Update solid particle velocity and position

for p = 1:spCount
    for i = 1:3
     xpid = CONNECTX{p}(i); % global x coordinate of basis function
    for j = 1:3
         ypid = CONNECTY{p}(j); % global y coordinate of basis function
         
         if nmass{i,j} ==0
             continue
         end
         
         v_ssp(p,:)                      = v_ssp(p,:) + dt * Nx{p}(i) * Ny{p}(j) * nforce{xpid,ypid}/nmass{xpid,ypid};
         x_sp(p,:)                       = x_sp(p,:) + dt * Nx{p}(i) * Ny{p}(j) * nmomentum{xpid,ypid}/nmass{xpid,ypid};
         d_sp(p,:)                       = x_sp(p,:) - x_spo(p,:);
    end
    end
end
% velocity particle: v_ssp
% position particle: x_sp
% displacement particle: d_sp

%% Mapping nodal velocity back to node
% Node variables
 nvelo                = cell(Nelemsx,Nelemsy);                   % Nodal Velocity
 nmomentum            = cell(Nelemsx,Nelemsy);                   % Nodal Momentum
 
 for i = 1:Nelemsx
    for j = 1:Nelemsy
        nmomentum{i,j}            = zeros(1,2);
        nvelo{i,j}              = zeros(1,2);  
    end
 end

 % Interpolation momentum from particle to node
 for p=1:spCount 
 for i = 1:3
     xpid = CONNECTX{p}(i); % global x coordinate of basis function
     for j = 1:3
         ypid = CONNECTY{p}(j); % global y coordinate of basis function
        
         % Momentum
         nmomentum{xpid,ypid}   = nmomentum{xpid,ypid} + m_sp(p) * v_ssp(p,:) * Nx{p}(i) * Ny{p}(j);
     end
 end
 end

 Momentum = 0;
 for i=1:Nelemsx
     for j = 1:Nelemsy
         Momentum = Momentum + nmomentum{i,j};
     end
 end
 
 % Velocity
 for i = 1:Nelemsx
    for j = 1:Nelemsy
        nvelo{i,j}     	= nmomentum{i,j} / nmass{i,j};
        
        % Boundary condition
        if i == 1 || i == Nelemsx
            nvelo{i,j}(1)      = 0;
        end
        
        if j == 1 || j == Nelemsy
            nvelo{i,j}(2)      = 0;
        end
    end
 end

 test = zeros(Nelemsx,1);
 for i = 1:Nelemsx
     test(i) = nvelo{i,2}(1);
 end
 
 
%% Update effective stress
% Calculate stress for solid phase
L_sp = cell(spCount,1);

for sp = 1:spCount
    L_sp{sp} = zeros(2,2);
end


for p = 1:spCount
    for i = 1:3
     xpid = CONNECTX{p}(i); % global x coordinate of basis function
    for j = 1:3
         ypid = CONNECTY{p}(j); % global y coordinate of basis function
         
         L_sp{p}   = L_sp{p} + nvelo{xpid,ypid}' * [dNx{p}(i)*Ny{p}(j) Nx{p}(i)*dNy{p}(j)];
    end
    end
    dESP = (L_sp{p} + L_sp{p}')/2*dt; 
        
        F_sp{p} = (eye(2,2)+L_sp{p}*dt)*F_sp{p};                           
        J = det(F_sp{p});
        V_sp(p)=V_spo(p)*J;   

        switch CModel
            case 'Neo_Hookean_Elastic'
                [s_sp(p,:)]=Neo_Hookean_elastic(CModel_parameter,F_sp{p},J);
            case 'Linear_Elastic'
                [s_sp(p,:)]=Linear_elastic(CModel_parameter,dESP,s_sp(p,:));             
            case 'Water'
                [s_sp(p,:)]=Water(CModel_parameter,(L_sp{p} + L_sp{p}')/2,J);
        end
        p_sp(p) = m_sp(p)/V_sp(p);
end      

test1 = 1;