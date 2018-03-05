function [nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si]=Interpolate_Particle_To_Grid(CONNECT,le,N,dG,spCount,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si)

%% Interpolation from particle to grid task
% Output
% nmass_si: nodal mass
% nmomentum_si: nodal momentum
% niforce_si: nodal internal force
% neforce_si: nodal external force
% traction_si: nodal traction

for sp=1:spCount
 % Build stress tensor
 SSP = [s_sp(sp,1) s_sp(sp,3);s_sp(sp,3) s_sp(sp,2)];
 
 for j=1:4
     npid                           = CONNECT(sp,j);
     
 % Mass
 nmass_si(npid)            = nmass_si(npid) + m_sp(sp)*N(sp,j);
 
 % Momentum
 nmomentum_si(npid,:)      = nmomentum_si(npid,:) + m_sp(sp)*v_ssp(sp,:)*N(sp,j);
 
 % Internal forces
niforce_si(npid,:)         = niforce_si(npid,:) - (V_sp(sp)*SSP*dG((sp-1)*4+j,:)')';

 % External forces
neforce_si(npid,:)         = neforce_si(npid,:) + b_sp*m_sp(sp)*N(sp,j);

% Traction
 traction_si(npid,:)       = traction_si(npid,:) + V_sp(sp)*ptraction_sp(sp,:)*N(sp,j)/le(1,1)/le(1,2);
end 
 end