function initial_figure = Plot_Initial(x_sp,LOC,le)

initial_figure=figure;
plot(x_sp(:,1),x_sp(:,2),'o');
title('initial condition');
grid on
axis([0,max(LOC(:,1)),0,max(LOC(:,2))]);
set(gca,'xtick',[0:le(1):max(LOC(:,1))]);
set(gca,'ytick',[0:le(2):max(LOC(:,2))]);
% set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);