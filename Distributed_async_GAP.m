%Async communication comparisons
clc
clear
close all
%Random Seed

%rng(3)%with c = randn(robotNum*taskNum,1); %Can lead to infeasible configurations
rng(10)
%Assignment Problem only ineq assignment constraints
%Num Robots/Agents
robotNum =50;
%Num Tasks
taskNum = robotNum-10;

%Relaxation constant close to one
xi = 1-10^-6;
%Random Objective function

% c = randn(robotNum*taskNum,1).^2;
c = rand(robotNum*taskNum,1);
%Fixed krocker formulation of constraints
A_assignment = kron(ones(robotNum,1),eye(taskNum,taskNum))';
A_local = kron(eye(robotNum,robotNum),ones(1,taskNum));
beq = ones(taskNum,1);
minLC = 0;
lb_relax = (0.5-xi)*ones(size(c,1),1);
ub_relax = 1+(xi-0.5)*ones(size(c,1),1);
lb = zeros(1,robotNum*taskNum);
ub = ones(1,robotNum*taskNum);
LocalConstraints =((minLC*abs(rand(size(A_local,2),1))+1).*A_local')';
Aineq_local = LocalConstraints;
maxLC = 1.2;
maxLocalConstraints = maxLC*ones(robotNum,1);
bineq = maxLocalConstraints;
%Enlarged Inner Parallel Set Computation
EIPS = -sum(abs([A_assignment;Aineq_local])')'./2+xi;

%Algorithm parameters & initialization
xInitial = ones(size(c,2),1);
lamInit = zeros(size([beq;bineq],1),1);
stepSize = [10^-3,10^-2];
tolerance = 10^-4;
maxIterations = 10^5;
ubPrimal =1+(xi-0.5);
lbPrimal =(0.5-xi);
lbDual = 0;
ubDual = 10;

%Compute Optimal solution using SCIP mixed integer solver
intcon = 1:robotNum*taskNum;
opts = optiset('solver','scip');

SCIPOptimization = opti('f',c,'ineq',[A_assignment;Aineq_local],[beq;bineq],'bounds',lb,ub,'xtype',intcon,'options',opts);
[xIntegerOptimal, fval_only_geq, exitflag, info] = solve(SCIPOptimization);


%Check if the relaxed problem is granular
opts = optiset('solver','scip');
SCIPOptimization = opti('f',zeros(size(c)),'ineq',[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS,'bounds',lb_relax,ub_relax,'options',opts);
[isGranularCheck, fval_granCheck, granular, info] = solve(SCIPOptimization);
%-1 means infeasible
%1 means feasible
if granular == 1
    disp("Granular");
else
    disp("Not Granular");
end

%Compute Regularization error terms
%Original
primalReg =10^-3; 
dualReg = 10^-1;

slaterTerms = (c'*isGranularCheck+(primalReg/2)*norm(isGranularCheck)^(2))/(-A_assignment(1,:)*isGranularCheck-beq(1)-EIPS(1));
regCorrectionTest = norm(A_assignment(1,:)*isGranularCheck-beq(1)-EIPS(1))*slaterTerms*sqrt(dualReg/(2*primalReg));


%Async Primal Dual Algorithm Parameters
scBlock = false;
primalNum = robotNum;
dualNum = robotNum;
commRate = 0.1;
plotVals = false;

%Async Primal Dual Algorithm
[x_async, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS+regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,plotVals,xIntegerOptimal);
% roundedGran_async = round(x_async);
diff_async = norm(x_async-xIntegerOptimal,1)/(1e-10+norm(xIntegerOptimal,1));

commRate = 0.5;
[x_async_0_5, es_async_0_5,numDualUpdates_async_0_5,convDist_async_0_5,constr_async_0_5,x_async_comparedtoOptimalOverTime_0_5] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS+regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,plotVals,xIntegerOptimal);
commRate = 0.75;
[x_async_0_75, es_async_0_75,numDualUpdates_async_0_75,convDist_async_0_75,constr_async_0_75,x_async_comparedtoOptimalOverTime_0_75] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS+regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,plotVals,xIntegerOptimal);
if all(Aineq_local*x_async<=bineq) && all(A_assignment*x_async<=beq)
    disp("Async Form is Feasible comm 0.1");
end
if all(Aineq_local*x_async_0_75<=bineq) && all(A_assignment*x_async_0_75<=beq)
    disp("Async Form is Feasible comm 0.75");
end
if all(Aineq_local*x_async_0_5<=bineq) && all(A_assignment*x_async_0_5<=beq)
    disp("Async Form is Feasible comm 0.5");
end


%Solve Regularized Lagrangian in a centralized algorithm
[x_centralized, lamGradTest,xComparedToOptimalOverTime,lamList] = centralized_PD_Alg(xInitial,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS+regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);
% roundedGran_centralized = round(x_centralized);

if all(Aineq_local*x_centralized<=bineq) && all(A_assignment*x_centralized<=beq)
    disp("Centralized Form is Feasible");
end



diff_centralized = norm(x_centralized-xIntegerOptimal,1)/(1e-10+norm(xIntegerOptimal,1));
diff_async2_central = norm(x_centralized-x_async);



%%
%Suboptimality Error Bound
% Hoffman(Aconst)
h = 1;
% kap = sqrt(p*q);
kap = 1;
apriori_optimalityEB = (norm(c,2)*(Hoffman([A_assignment;Aineq_local])*norm([A_assignment;Aineq_local],inf)+kap)+norm(c,2)*kap)*(h/2);
apriori_feasibilityEB = (Hoffman([A_assignment;Aineq_local])*norm([A_assignment;Aineq_local],inf)+2*kap)*(h/2);

%A posteriori assumes we know the optimal solution centralized
%Note: I am taking the 2-norm for distance instead of the infimum
apost_optimalityEB_central = norm(c,2)*(norm(x_centralized-xIntegerOptimal,2)+kap*0.5)+norm(c,2)*kap*(h/2);
apost_feasibilityEB_central = norm(x_centralized-xIntegerOptimal,2)+kap*h;

%Async form
apost_optimalityEB_async = norm(c,2)*(norm(x_async-xIntegerOptimal,2)+kap*0.5)+norm(c,2)*kap*(h/2);
apost_feasibilityEB_async = norm(x_async-xIntegerOptimal,2)+kap*h;
%%

%Plotting differences in communication rates
%  figure()
% 
% hold on
% semilogy(es1(2:500)','-','LineWidth',2)
% semilogy(es075(2:500)','-','LineWidth',2)
% semilogy(es05(2:500)','-','LineWidth',2)
% semilogy(es01(2:500)','-','LineWidth',2)
% % semilogy(esLP(2:end)','-','LineWidth',2)
% title('Cost')
% xlabel('Iteration Number','FontWeight','Bold')
% ylabel('$c^{T}\cdot z(k)$','Interpreter','latex','FontWeight','Bold')
% hold off
% legend('Comm Rate=1.0','Comm Rate=0.75','Comm Rate=0.5','Comm Rate=0.1');
% %%
% figure ()
% fontSizeOverall = 25;
% hold on
% semilogy(convDist1(2:300)','-','LineWidth',2)
% semilogy(convDist075(2:300)','-','LineWidth',2)
% semilogy(convDist05(2:300)','-','LineWidth',2)
% semilogy(convDist01(2:300)','-','LineWidth',2)
% % semilogy(convDistLP(2:500)','-','LineWidth',2)
% ax = gca;
% ax.FontSize = fontSizeOverall;
% % semilogy(convDistLP(2:end)','-','LineWidth',2)
% title('Distance Between Iterates Convergence Comparison','FontSize',fontSizeOverall)
% xlabel('Iteration Number','FontWeight','Bold','FontSize',fontSizeOverall)
% ylabel('$||z(k) - z(k-1)||$','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
% hold off
% legend('Comm Rate=1.0','Comm Rate=0.75','Comm Rate=0.5','Comm Rate=0.1','FontSize',fontSizeOverall);
% 
% %%
% 
% 
% 
% figure()
% % 
% hold on
% semilogy(constr1(2:500)','-','LineWidth',2)
% semilogy(constr075(2:500)','-','LineWidth',2)
% semilogy(constr05(2:500)','-','LineWidth',2)
% semilogy(constr01(2:500)','-','LineWidth',2)
% semilogy(constrLP(2:500)','-','LineWidth',2)
% title('Constraint Satisfaction Relative to Iteration Number')
% xlabel('Iteration Number','FontWeight','Bold')
% ylabel('$Az(k)-b$','Interpreter','latex','FontWeight','Bold')
% hold off
% legend('Comm Rate=1.0','Comm Rate=0.75','Comm Rate=0.5','Comm Rate=0.1');



%%
close all
figure ()
fontSizeOverall = 25;
hold on
plotLen = 5000;
if length(xComparedToOptimalOverTime)>= plotLen
    semilogy(xComparedToOptimalOverTime(2:plotLen)','-','LineWidth',2)
    semilogy(x_async_comparedtoOptimalOverTime_0_75(2:plotLen)','-','LineWidth',2)
    semilogy(x_async_comparedtoOptimalOverTime_0_5(2:plotLen)','-','LineWidth',2)
    semilogy(x_async_comparedtoOptimalOverTime(2:plotLen)','-','LineWidth',2)

else
    semilogy(xComparedToOptimalOverTime(2:end)','-','LineWidth',2)
    semilogy(x_async_comparedtoOptimalOverTime_0_75(2:end)','-','LineWidth',2)
    semilogy(x_async_comparedtoOptimalOverTime_0_5(2:end)','-','LineWidth',2)
    semilogy(x_async_comparedtoOptimalOverTime(2:end)','-','LineWidth',2)
end

%Plot Theoretical Error Bounds
EB_plot_length = length(x_async_comparedtoOptimalOverTime(2:end));
eb_xPlot = linspace(0,EB_plot_length,EB_plot_length);
eb_yPlot = apost_feasibilityEB_async *ones(EB_plot_length);
plot(eb_xPlot,eb_yPlot,'g')
title('Distance to Optimum','FontSize',fontSizeOverall)
ylabel('$||z_{\kappa}(k) - z^{*}_{MILP}||$','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
xlabel('Iteration Number','FontWeight','Bold','FontSize',fontSizeOverall)
ax = gca;
ax.FontSize = fontSizeOverall;
% ylabel('Distance Between Iterates','FontWeight','Bold')
% hold off
legend('Comm Rate=1.0','Comm Rate=0.75','Comm Rate=0.5','Comm Rate=0.1','A-Post Feasibility','Interpreter','latex','FontSize',fontSizeOverall);
% legend('Comm Rate=1.0','Comm Rate=0.1','Interpreter','latex','FontSize',fontSizeOverall);

