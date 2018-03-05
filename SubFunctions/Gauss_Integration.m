function [nmass_si,nmomentum_si,niforce_si] = Gauss_Integration(Reconstruction_density,Reconstruction_momentumx,Reconstruction_momentumy,Reconstruction_stressx,Reconstruction_stressy,Reconstruction_shear,CONNECT,LOC,le,GaussCount,x_gauss,y_gauss,omega_gauss,spElems,active_elements,nmass_si,nmomentum_si,niforce_si)

% Compute the Gauss Integration for each elements

ElemsCount = length(active_elements);
for i=1:ElemsCount
    k       = active_elements(i);
    spid    = find(spElems==k);
    for g=1:GaussCount
        position_Gauss = [x_gauss(i,g);y_gauss(i,g)];
        for j=1:4
            npid    = CONNECT(spid(1),j);
            [N_gauss,dNx_gauss,dNy_gauss]=linearshape(position_Gauss,LOC(npid,:),le(1),le(2));
            % Mass
            nmass_si(npid)            = nmass_si(npid) + Reconstruction_density(i,g)*omega_gauss(i,g)*N_gauss; 
            % Momentum
            nmomentum_si(npid,1)      = nmomentum_si(npid,1) + Reconstruction_momentumx(i,g)*omega_gauss(i,g)*N_gauss; 
            nmomentum_si(npid,2)      = nmomentum_si(npid,2) + Reconstruction_momentumy(i,g)*omega_gauss(i,g)*N_gauss;
            % Stress
            % Build stress tensor
            SSP = [Reconstruction_stressx(i,g) Reconstruction_shear(i,g);Reconstruction_shear(i,g) Reconstruction_stressy(i,g)];
            niforce_si(npid,:)        = niforce_si(npid,:)  - [dNx_gauss dNy_gauss] * SSP * omega_gauss(i,g);
        end
    end 
end