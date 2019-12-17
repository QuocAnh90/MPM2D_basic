function [nvelo_si] = Interpolate_velocity_back_MLS(NODES,nodeCount,spCount,CONNECT,v_ssp,N_MLS)

%% Interpolation from particle to grid task
% Node variables
nvelo_si                = zeros(nodeCount,2);                   % Nodal Velocity

for p = 1 : spCount

     for j=1 : NODES(p)
         npid                           = CONNECT{p}(j);

              if N_MLS{p}(j)==0
             continue
              end

     % Velo
     nvelo_si(npid,:)          = nvelo_si(npid,:) + v_ssp(p,:)*N_MLS{p}(j);     
     end
end