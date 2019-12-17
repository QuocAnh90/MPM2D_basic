function [s_sp]=Water(CModel_parameter,dESP,J)

% Input
% K: bulk modulus
% u: vicosity
% gamma: rate
K = CModel_parameter(1);
u = CModel_parameter(2);
gamma = CModel_parameter(3);


% Rate of deformation tensor
D = dESP;

% The Deviatoric part of the rate of deformation tensor
Dprime = D - 1/3*trace(D)*eye(2,2);
Shear = 2*u*Dprime;

% The isotropic pressure
p = K*(J^gamma-1);

% Cauchy stress 
SSP = p*eye(2,2) + Shear;

% Regroup the stress
s_sp(1) = SSP(1,1);
s_sp(2) = SSP(2,2);
s_sp(3) = SSP(1,2);

