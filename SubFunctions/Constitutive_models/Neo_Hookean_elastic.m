function [s_sp]=Neo_Hookean_elastic(CModel_parameter,F_sp,J)

% This is the hyper-elastic Neo-Hookean model

% Parameter of the model
E = CModel_parameter(1);
nu = CModel_parameter(2);

Lambda = E*nu/(1+nu)/(1-2*nu);
Mu     = E/2/(1+nu);

SSP = Lambda*log(J)/J*eye(2,2) + Mu/J*(F_sp*F_sp'-eye(2,2));
        s_sp(1) = SSP(1,1);
        s_sp(2) = SSP(2,2);
        s_sp(3) = SSP(1,2);
        