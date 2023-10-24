%Async communication comparisons
clc
clear
close all
%Random Seed

%rng(3)%with c = randn(robotNum*taskNum,1); %Can lead to infeasible configurations
rng(7)
%Assignment Problem only ineq assignment constraints
%Num Robots/Agents
robotNum =100;
%Num Tasks
taskNum = robotNum;

%Relaxation constant close to one
xi = 1-10^-3;
%Random Objective function

c = randn(robotNum*taskNum,1);

%Fixed krocker formulation of constraints
A_assignment = kron(ones(robotNum,1),eye(taskNum,taskNum))';
A_local = kron(eye(robotNum,robotNum),ones(1,taskNum));
beq = ones(taskNum,1);
minLC = 0;
lb_relax = (0.5-xi)*ones(size(c,1),1);
ub_relax = 1+(xi-0.5)*ones(size(c,1),1);
lb = zeros(1,robotNum*taskNum);
ub = ones(1,robotNum*taskNum);
LocalConstraints =((minLC*abs(round(rand(size(A_local,2),1))+1)).*A_local')';
Aineq_local = LocalConstraints;
maxLC = 2.2;
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
primalReg =10^-5; 
dualReg = 10^-4;
% primalReg =10^-2; 
% dualReg = 10^-1;


slaterTerms = (c'*isGranularCheck+(primalReg/2)*norm(isGranularCheck)^(2))/(-A_assignment(1,:)*isGranularCheck-beq(1)-EIPS(1));
regCorrectionTest = norm(A_assignment(1,:)*isGranularCheck-beq(1)-EIPS(1))*slaterTerms*sqrt(dualReg/(2*primalReg));

%%
% Testing Slater's
% [x_centralized, lamGradTest,xComparedToOptimalOverTime,lamList] = centralized_PD_Alg(xInitial,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);
[x_centralized, lamGradTest,xComparedToOptimalOverTime,lamList] = centralized_PD_Alg(xInitial,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);
% roundedGran_centralized = round(x_centralized);

if all(Aineq_local*x_centralized<=bineq) && all(A_assignment*x_centralized<=beq)
    disp("Centralized Form is Feasible");
end


%%
% Fixing minor infeasibility with a bit of low async
commRate = 1.0;
computeRate = 0.5;
maxIterations = 10^2;
[x_async_075, ~,~,~,~,~] = Async_PD(x_async_075,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
if sum(Aineq_local*x_async_0_5<=bineq) && all(A_assignment*x_async_0_5<=beq)
    disp("Original Comp of 0.75 is now feasible after being used to warm start Comm of 0.5");
end
% [x_async_0_5, ~,~,~] = centralized_PD_Alg(x_async_0_5,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);
%%
[x_async_0_1, ~,~,~] = centralized_PD_Alg(x_async_0_1,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);

%%
%Fixing minor infeasibility with a bit of low async
[x_async_05, ~,~,~] = centralized_PD_Alg(x_async_05,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);
%%
%Async Primal Dual Algorithm Parameters
scBlock = false;
primalNum = robotNum;
dualNum = robotNum;
commRate = 0.1;
plotVals = false;
compute_rate = 1.0;



%Async Primal Dual Algorithm
% Testing if different vals are needed
maxIterations = 10^4;

commRate = 0.5;
[x_async_0_5, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
% roundedGran_async = round(x_async);
diff_async = norm(x_async_0_5-xIntegerOptimal,1)/(1e-10+norm(xIntegerOptimal,1));

if all(Aineq_local*x_async_0_5<=bineq) && all(A_assignment*x_async_0_5<=beq)
    disp("Async Form is Feasible comm 0.5");
end

%%
commRate = 0.1;
[x_async, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
commRate = 0.5;
[x_async_0_5, es_async_0_5,numDualUpdates_async_0_5,convDist_async_0_5,constr_async_0_5,x_async_comparedtoOptimalOverTime_0_5] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
commRate = 0.75;
[x_async_0_75, es_async_0_75,numDualUpdates_async_0_75,convDist_async_0_75,constr_async_0_75,x_async_comparedtoOptimalOverTime_0_75] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
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
[x_centralized, lamGradTest,xComparedToOptimalOverTime,lamList] = centralized_PD_Alg(xInitial,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);
% roundedGran_centralized = round(x_centralized);

if all(Aineq_local*x_centralized<=bineq) && all(A_assignment*x_centralized<=beq)
    disp("Centralized Form is Feasible");
end

%%

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
close all
figure ()
fontSizeOverall = 25;
hold on
plotLen = 1000000;
% if length(xComparedToOptimalOverTime)>= plotLen 
%     semilogy(xComparedToOptimalOverTime(2:plotLen)','-','LineWidth',2)
%     semilogy(x_async_comparedtoOptimalOverTime_0_75(2:plotLen)','-','LineWidth',2)
%     semilogy(x_async_comparedtoOptimalOverTime_0_5(2:plotLen)','-','LineWidth',2)
%     semilogy(x_async_comparedtoOptimalOverTime(2:plotLen)','-','LineWidth',2)
% 
% else
    semilogy(xComparedToOptimalOverTime(2:end)','-','LineWidth',2)
    semilogy(x_async_comparedtoOptimalOverTime_0_75(2:end)','-','LineWidth',2)
    semilogy(x_async_comparedtoOptimalOverTime_0_5(2:end)','-','LineWidth',2)
    semilogy(x_async_comparedtoOptimalOverTime(2:end)','-','LineWidth',2)
% end

%Plot Theoretical Error Bounds
% EB_plot_length = length(x_async_comparedtoOptimalOverTime(2:plotLen));
% eb_xPlot = linspace(0,EB_plot_length,EB_plot_length);
% eb_yPlot = apost_feasibilityEB_async *ones(EB_plot_length);
% plot(eb_xPlot,eb_yPlot,'g')
title('Distance to Optimum','FontSize',fontSizeOverall)
ylabel('$||z_{\kappa}(k) - z^{*}_{MILP}||$','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
xlabel('Iteration Number','FontWeight','Bold','FontSize',fontSizeOverall)
ax = gca;
ax.FontSize = fontSizeOverall;
% ylabel('Distance Between Iterates','FontWeight','Bold')
% hold off
legend('Comm Rate=1.0','Comm Rate=0.75','Comm Rate=0.5','Comm Rate=0.1','A-Post Feasibility','Interpreter','latex','FontSize',fontSizeOverall);
% legend('Comm Rate=1.0','Comm Rate=0.1','Interpreter','latex','FontSize',fontSizeOverall);

%%
%Testing varying computation rate
fontSizeOverall = 25;
% close all
%Async Primal Dual Algorithm Parameters
scBlock = false;
primalNum = robotNum;
dualNum = robotNum-30;
commRate = 1.0;
plotVals = false;
compute_rate = 1.0;
maxIterations = 10^5;
tolerance = 10^-5;
%Async Primal Dual Algorithm
%spmd
[x_async, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
all(Aineq_local*x_async<=bineq) && all(A_assignment*x_async<=beq)
compute_rate = 0.75;
[x_async_075, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime_075] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);

all(Aineq_local*x_async_075<=bineq) && all(A_assignment*x_async_075<=beq)
compute_rate = 0.5;
[x_async_05, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime_05] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
all(Aineq_local*x_async_05<=bineq) && all(A_assignment*x_async_05<=beq)
compute_rate = 0.1;
[x_async_01, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime_01] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
all(Aineq_local*x_async_01<=bineq) && all(A_assignment*x_async_01<=beq)

[x_centralized, lamGradTest,xComparedToOptimalOverTime,lamList] = centralized_PD_Alg(xInitial,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);
%end
%%
fontSizeOverall = 25;
figure()
hold on
semilogy(x_async_comparedtoOptimalOverTime(2:end)','-','LineWidth',2)
semilogy(x_async_comparedtoOptimalOverTime_075(2:end)','-','LineWidth',2)
semilogy(x_async_comparedtoOptimalOverTime_05(2:end)','-','LineWidth',2)
semilogy(x_async_comparedtoOptimalOverTime_01(2:end)','-','LineWidth',2)
%Plot Theoretical Error Bounds
% EB_plot_length = length(x_async_comparedtoOptimalOverTime(2:plotLen));
% eb_xPlot = linspace(0,EB_plot_length,EB_plot_length);
% eb_yPlot = apost_feasibilityEB_async *ones(EB_plot_length);
% plot(eb_xPlot,eb_yPlot,'g')
ax = gca;
ax.FontSize = fontSizeOverall;
titleFull = strcat('Distance to Optimum');
title(titleFull,'FontSize',fontSizeOverall)
ylabel('$||z_{\kappa}(k) - z^{*}_{MILP}||$','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)

xlabel('Iteration Number','FontWeight','Bold','FontSize',fontSizeOverall)
legend('Comp Rate=1.0','Comp Rate=0.75','Comp Rate=0.5','Comp Rate=0.1','A-Post Feasibility','Interpreter','latex','FontSize',fontSizeOverall);
ylim([0, inf])

%%
%Finding regularization terms that work

primalReg =10^-5; 
dualReg = 10^-4;
% primalReg =10^-2; 
% dualReg = 10^-1;

slaterTerms = (c'*isGranularCheck+(primalReg/2)*norm(isGranularCheck)^(2))/(-A_assignment(1,:)*isGranularCheck-beq(1)-EIPS(1));
regCorrectionTest = norm(A_assignment(1,:)*isGranularCheck-beq(1)-EIPS(1))*slaterTerms*sqrt(dualReg/(2*primalReg));

opts = optiset('solver','scip');
SCIPOptimization = opti('f',c,'ineq',[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,'bounds',lb_relax,ub_relax,'options',opts);
[testing123, fval_granCheck, granular, info] = solve(SCIPOptimization);
scipTest = round(testing123);

all(Aineq_local*scipTest<=bineq) && all(A_assignment*scipTest<=beq)
