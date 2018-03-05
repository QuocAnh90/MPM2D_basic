function [x_gauss,y_gauss,omega_gauss,active_elements] = Locate_Gauss_Point(GaussCount,spElems,LOCC,le)

%% copy the code written by Lisa Wobbes, 27 Feb 2018
% Compute the Gauss point of the element with weights

% 4 Gauss Points recommended
% Gauss Point in the element: 1=bottomleft, 2=bottomright, 3= topright, 4=top left
active_elements = unique(spElems);
ElemsCount = length (active_elements);

x_gauss = zeros(ElemsCount, GaussCount);
y_gauss = zeros(ElemsCount, GaussCount);
omega_gauss = zeros(ElemsCount, GaussCount);

for i=1:ElemsCount
    k = active_elements(i);                 % Index of the element
    
    x0 = LOCC(k,1)-0.5*le(1);               
    x1 = LOCC(k,1)+0.5*le(1);   
    [x_gauss_temp, x_omega_gauss] = lgwt(GaussCount/2,x0,x1);
    x_gauss(i,:) = [x_gauss_temp(1), x_gauss_temp(2), x_gauss_temp(2), x_gauss_temp(1)];
    
    y0 = LOCC(k,2)-0.5*le(2);               
    y1 = LOCC(k,2)+0.5*le(2);   
    [y_gauss_temp, y_omega_gauss] = lgwt(GaussCount/2,y0,y1);
    y_gauss(i,:) = [y_gauss_temp(1), y_gauss_temp(1), y_gauss_temp(2), y_gauss_temp(2)];
    
    omega_gauss_temp = x_omega_gauss*y_omega_gauss';
    omega_gauss(i,:) = reshape(omega_gauss_temp, 1, 4);
end