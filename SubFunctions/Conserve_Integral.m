function [Conserve_mass,Conserve_momentum] = Conserve_Integral(le,active_elements,mspoints,m_sp,v_ssp)

%% copy the code written by Lisa Wobbes, 27 Feb 2018
% Compute the f

ElemsCount = length(active_elements);
element_size = le(1)*le(2);

Conserve_mass = zeros(ElemsCount,1);
Conserve_momentum = zeros(ElemsCount,2);

for i=1:ElemsCount
    k = active_elements(i);                 % Index of the element
    for p=1:length(mspoints{k})
        pid = mspoints{k}(p);
        Conserve_mass(i) = Conserve_mass(i) + m_sp(pid);
        Conserve_momentum(i,:) = Conserve_momentum(i,:) + m_sp(pid)*v_ssp(pid,:);
    end
    Conserve_mass(i) = Conserve_mass(i)/element_size;
    Conserve_momentum(i,:) = Conserve_momentum(i,:)/element_size;
end

