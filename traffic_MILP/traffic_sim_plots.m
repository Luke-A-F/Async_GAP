%Note exportgraphics as pdf with the following command
%exportgraphics(gcf, 'soln_limit_dynamics.pdf');


clc
clear all
close all
load("traffic_MILP_comm_1_0.mat")
load("traffic_MILP_comm_075.mat")
load("traffic_MILP_comm_05.mat")
load("traffic_MILP_comm_01.mat")



%%
fontSizeOverall = 25;
axis.FontSize =fontSizeOverall;

hold on


semilogy(constr_async_1_0,"LineWidth", 2);
semilogy(constr_async_075,"LineWidth", 2);
semilogy(constr_async_05,"LineWidth", 2);
semilogy(constr_async_01,"LineWidth", 2);

ax = gca;
% xlim([0 10^4])
title("Feasibility")
ylabel("Constraints Violated",'Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
xlabel("Iteration Number",'Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
legend("Comm Rate=1.0","Comm Rate=0.75","Comm Rate=0.5","Comm Rate=0.1");
ax.FontSize = fontSizeOverall;

figure()
hold on
x_async_comparedtoOptimalOverTime_1_0 = x_async_comparedtoOptimalOverTime_1_0';
x_async_comparedtoOptimalOverTime_075 = x_async_comparedtoOptimalOverTime_075';
x_async_comparedtoOptimalOverTime_05 = x_async_comparedtoOptimalOverTime_05';
x_async_comparedtoOptimalOverTime_01 = x_async_comparedtoOptimalOverTime_01';

semilogy(x_async_comparedtoOptimalOverTime_1_0,"LineWidth", 2);
semilogy(x_async_comparedtoOptimalOverTime_075,"LineWidth", 2);
semilogy(x_async_comparedtoOptimalOverTime_05,"LineWidth", 2);
semilogy(x_async_comparedtoOptimalOverTime_01,"LineWidth", 2);

ax = gca;
% xlim([0 10^4])

title("Distance to Optimum")
ylabel({'$||z_{\kappa}(k) - z^{*}_{MILP}||$'},'Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
% ylabel('||z(k) - z^{*}_{MILP}||','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
xlabel("Iteration Number",'Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
legend("Comm Rate=1.0","Comm Rate=0.75","Comm Rate=0.5","Comm Rate=0.1");
ax.FontSize = fontSizeOverall;

figure()
hold on


semilogy(es_async_1_0,"LineWidth", 2);
semilogy(es_async_075,"LineWidth", 2);
semilogy(es_async_05,"LineWidth", 2);
semilogy(es_async_01,"LineWidth", 2);

ax = gca;
% xlim([0 10^4])

title("Relative Optimality Gap")
ylabel({'$\frac{|c^{T}z^{*}_{MILP}- c^{T}[z_{\kappa}(k)]_{r}|}{|c^{T}z^{*}_{MILP}|}$'},'Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
% ylabel('||z(k) - z^{*}_{MILP}||','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
xlabel("Iteration Number",'Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)

legend("Comm Rate=1.0","Comm Rate=0.75","Comm Rate=0.5","Comm Rate=0.1");
ax.FontSize = fontSizeOverall;

