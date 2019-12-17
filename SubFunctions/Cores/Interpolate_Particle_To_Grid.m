function [nmass,nmomentum,niforce,neforce,traction]=Interpolate_Particle_To_Grid(NODES,nodeCount,CONNECT,le,N,dN,spCount,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m)

%% Interpolation from particle to grid task
% Node variables
nmass                = zeros(nodeCount,1);                   % Nodal Mass
nmomentum            = zeros(nodeCount,2);                   % Nodal Momentum
niforce              = zeros(nodeCount,2);                   % Nodal Internal force
neforce              = zeros(nodeCount,2);                   % Nodal External force
traction             = zeros(nodeCount,2);                   % Nodal Traction

 for p=1:spCount
 % Build stress tensor
 SSP = [s_sp(p,1) s_sp(p,3);s_sp(p,3) s_sp(p,2)];
 
 for j=1:NODES(p)
     npid                           = CONNECT{p}(j);
     
          if N{p}(j)==1
         continue
          end
     
 % Mass
 nmass(npid)            = nmass(npid) + m(p)*N{p}(j);
 
 % Momentum
 nmomentum(npid,:)      = nmomentum(npid,:) + m(p)*v_ssp(p,:)*N{p}(j);
 
 % Internal forces
niforce(npid,:)         = niforce(npid,:) - (V_sp(p)*SSP*dN{p}(:,j))';

 % External forces
neforce(npid,:)         = neforce(npid,:) + b_sp(p,:)*m(p)*N{p}(j);

% Traction
%  traction_si(npid,:)       = traction_si(npid,:) + V_sp(sp)*ptraction_sp(sp,:)*N{sp}(j)/le(1,1)/le(1,2);
 traction(npid,1)       = traction(npid,1) + V_sp(p)*ptraction_sp(p,1)*N{p}(j)/le(1);
 traction(npid,2)       = traction(npid,2) + V_sp(p)*ptraction_sp(p,2)*N{p}(j)/le(2); 
 end 
 end
 
% Notes on traction calculation

% If traction / le: it means that the traction layer is equal to the cell
% layer
% If traction / lp: it means that the traction layer is equal to the
% particle  layer