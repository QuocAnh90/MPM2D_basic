function [nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si]=Interpolate_Particle_To_Grid_MLS(NODES,GAUSS,NODES_G,cellCount,nodeCount,...
    CONNECT,CONNECT_G,CONNECT_C,...
    le,N,N_g,N_c,dN_c,N_MLS,...
    spCount,b_sp,V_sp,p_sp,ptraction_sp,v_ssp,s_sp,m_sp)

%% Interpolation from particle to grid task
% Node variables
nmass_si                = zeros(nodeCount,1);                   % Nodal Mass
nvelo_si                = zeros(nodeCount,2);                   % Nodal Velocity
niforce_si              = zeros(nodeCount,2);                   % Nodal Internal force
neforce_si              = zeros(nodeCount,2);                   % Nodal External force
traction_si             = zeros(nodeCount,2);                   % Nodal Traction
p_sc                    = zeros(cellCount,1); 
SSC                     = cell(cellCount,1);
s_sc                    = zeros(cellCount,3);
b_sc                    = zeros(cellCount,2); 
Nt                      = zeros(cellCount,1);

% Cauchy tress tensor at cell
for c = 1:cellCount
    SSC{c} = zeros(2,2);
end

%% Mapping from particle to centroids
 for p = 1 : spCount
     % Build stress tensor
     SSP = [s_sp(p,1) s_sp(p,3);s_sp(p,3) s_sp(p,2)];

     % Mapping from particle to node
     for j=1 : NODES(p)
         npid                           = CONNECT{p}(j);

              if N{p}(j)==0
             continue
              end

     % Velo
     nvelo_si(npid,:)          = nvelo_si(npid,:) + v_ssp(p,:)*N_MLS{p}(j);     
     
     % Traction
     traction_si(npid,1)       = traction_si(npid,1) + V_sp(p)*ptraction_sp(p,1)*N{p}(j)/le(1);
     traction_si(npid,2)       = traction_si(npid,2) + V_sp(p)*ptraction_sp(p,2)*N{p}(j)/le(2); 
     end
     
     % Mapping from particle to centroids
     for c = 1:GAUSS(p)
     cpid = CONNECT_G{p}(c);
     
     % Density
     p_sc(cpid) = p_sc(cpid) + p_sp(p) * N_g{p}(c);
%      % Body
     b_sc(cpid,:) = b_sc(cpid,:) + p_sp(p) * b_sp(p,:) * N_g{p}(c);
     % Stress
     SSC{cpid} = SSC{cpid} + SSP * N_g{p}(c);
     
     Nt(cpid) = Nt(cpid) + N_g{p}(c);
    end
 end
 
 % Quadratic bspline mapping
 for c = 1:cellCount
     if Nt(c)==0
         continue
     end
      p_sc(c) = p_sc(c)/Nt(c);
      b_sc(c,:) = b_sc(c,:)/Nt(c);
      SSC{c}  = SSC{c}/Nt(c);
 end
 
 
 % Reorder stress tensor
 for c=1:cellCount
     s_sc(c,1) = SSC{c}(1,1);
     s_sc(c,2) = SSC{c}(2,2);
     s_sc(c,3) = SSC{c}(1,2);
 end
 
 %% Mapping from centroids to node
 for c = 1:cellCount
     if p_sc(c)==0
        continue
     end
     
     for j = 1:NODES_G(c);
         npid = CONNECT_C{c}(j);
             
     % Mass
     nmass_si(npid)            = nmass_si(npid) + p_sc(c) * N_c{c}(j) * le(1)*le(2);
     % Internal forces
     niforce_si(npid,:)         = niforce_si(npid,:) - (SSC{c} * dN_c{c}(:,j))' * le(1)*le(2);
     % External forces
     neforce_si(npid,:)         = neforce_si(npid,:) + b_sc(c,:) * N_c{c}(j) * le(1)*le(2);
     end
 end
 
 % Compute momentum
 nmomentum_si = zeros(nodeCount,2);
for n=1:nodeCount
    nmomentum_si(n,:) = nvelo_si(n,:) * nmass_si(n);
end

%% Testing momentum
momentum_p              = zeros(spCount,2);
for p=1:spCount
    momentum_p(p,:) = v_ssp(p,:) * m_sp(p);
end

 
 test=1;