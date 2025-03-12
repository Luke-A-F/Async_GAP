clc
clear all
fontSizeOverall = 25;
x = linspace(0,3);
y0 = zeros(length(x),1);
y1 = ones(length(x),1);
y2 = 2*ones(length(x),1);

%MILP Constraints
%y<=2.1-x
y_MILP = 2.1-x;
%xi = 0.8
%y<=2.3-x
y_xi_0_8 = 2.3-x;
%xi = 1
%y<=2.5-x
y_xi_1 = 2.5-x;

% region_MILP 
xlim([0,3]);
ylim([0,3]);
hold on
area(x,y_xi_1,'FaceColor','m','FaceAlpha',0.5,'EdgeAlpha',0.3,'LineWidth',2);
area(x,y_xi_0_8,'FaceColor','g','FaceAlpha',0.5,'EdgeAlpha',0.3,'LineWidth',2);
area(x,y_MILP,'FaceColor','b','FaceAlpha',0.5,'EdgeAlpha',0.3,'LineWidth',2);
% plot(x,y_MILP);

% plot(x,y_xi_0_8);

% plot(x,y_xi_1);
ax = gca;
ax.FontSize = fontSizeOverall;
plot(x,y0,'k','LineWidth',2);
plot(x,y1,'k','LineWidth',2);
plot(x,y2,'k','LineWidth',2);
xlabel('x','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall);
ylabel('y','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall);

legend('$\xi=1$','$\xi = 0.8$', 'MILP','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall);