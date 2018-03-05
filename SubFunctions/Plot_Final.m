function StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,e_sp,v_ssp,spCount)

% displacement = zeros(spCount,1);
velocity = zeros(spCount,1);
% strain = zeros(spCount,1);
for sp=1:spCount
%     displacement(sp) = sqrt(d_sp(sp,1)^2+d_sp(sp,2)^2);
    velocity(sp) = sqrt(v_ssp(sp,1)^2+v_ssp(sp,2)^2);
%     strain(sp) = sqrt(e_sp(sp,1)^2+e_sp(sp,2)^2);
end

StressProfile1=figure;
% set(StressProfile1, 'visible','off');
sz = 40;
color = velocity;
scatter(x_sp(:,1),x_sp(:,2),sz,color,'filled');
grid on
axis([0,max(LOC(:,1)),0,max(LOC(:,2))]);
set(gca,'xtick',[0:le(1):max(LOC(:,1))]);
set(gca,'ytick',[0:le(2):max(LOC(:,2))]);
h=colorbar;
colormap(jet(256))
% set(h, 'ytick', [0:0.2:1.2]);
    caxis([0 5]);