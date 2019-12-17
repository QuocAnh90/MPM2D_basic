function [s_sp]=Linear_elastic(CModel_parameter,dESP,s_sp)

% de_sp                   = zeros(spCount,3);                     % Strain increment
% e_sp                    = zeros(spCount,3);                     % Strain tensor
% ds_sp                   = zeros(spCount,3);                     % Stress increment


% Parameter of the model
E = CModel_parameter(1);
nu = CModel_parameter(2);

% Elastic matrix
D = [E/(1-nu.^2) E*nu/(1-nu.^2) 0; E*nu/(1-nu.^2) E/(1-nu.^2) 0; 0 0 E/(1-nu)];

% THis is the linear elastic model
        de_sp(1) = dESP(1,1);
        de_sp(2) = dESP(2,2);
        de_sp(3) = dESP(2,1);
        
        ds_sp = D * de_sp';
%         e_sp(spid,:)  =  e_sp(spid,:)+de_sp(spid,:);
        s_sp = s_sp + ds_sp';
        