%Async communication comparisons
clc
clear
close all
delete(gcp);
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
disp("SCIP Done");

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

slaterTerms = (c'*isGranularCheck+(primalReg/2)*norm(isGranularCheck)^(2))/(-A_assignment(1,:)*isGranularCheck-beq(1)-EIPS(1));
regCorrectionTest = norm(A_assignment(1,:)*isGranularCheck-beq(1)-EIPS(1))*slaterTerms*sqrt(dualReg/(2*primalReg));

%% Parallel Tests spmd tests
scBlock = false;
primalNum = robotNum;
dualNum = robotNum-30;
commRate = 0.1;
plotVals = false;
compute_rate = 1.0;
% p = parpool(4);
spmd(4)
%     labindex
% % spmdIndex = 1;
   [x_async,x_async_comparedtoOptimalOverTime] = parallelFunctionForSims(labindex,c,A_assignment,Aineq_local,beq,bineq,EIPS,regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
end
% delete(gcp);

x_async_1 = cell2mat(x_async(1,1));
x_async_0_75 = cell2mat(x_async(1,2));
x_async_0_5 = cell2mat(x_async(1,3));
x_async_0_1 = cell2mat(x_async(1,4));

xComparedToOptimalOverTime_1 = cell2mat(x_async_comparedtoOptimalOverTime(1,1));
x_async_comparedtoOptimalOverTime_0_75 = cell2mat(x_async_comparedtoOptimalOverTime(1,2));
x_async_comparedtoOptimalOverTime_0_5 = cell2mat(x_async_comparedtoOptimalOverTime(1,3));
x_async_comparedtoOptimalOverTime_0_1 = cell2mat(x_async_comparedtoOptimalOverTime(1,4));


save('parallelHundredvaryingCommRates')

% [x_centralized, lamGradTest,xComparedToOptimalOverTime,lamList] = centralized_PD_Alg(xInitial,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);
if all(Aineq_local*x_async_0_1<=bineq) && all(A_assignment*x_async_0_1<=beq)
    disp("Async Form is Feasible comm 0.1"); 
end
if all(Aineq_local*x_async_0_75<=bineq) && all(A_assignment*x_async_0_75<=beq)
    disp("Async Form is Feasible comm 0.75");
end
if all(Aineq_local*x_async_0_5<=bineq) && all(A_assignment*x_async_0_5<=beq)
    disp("Async Form is Feasible comm 0.5");
end
if all(Aineq_local*x_async_1<=bineq) && all(A_assignment*x_async_1<=beq)
    disp("Centralized Form is Feasible");
end
%% Parallel Tests spmd tests
%Varying computations
scBlock = false;
primalNum = robotNum;
dualNum = robotNum-30;
commRate = 1.0;
plotVals = false;
compute_rate = 1.0;
% p = parpool(4);
spmd(4)
%     labindex
% % spmdIndex = 1;
   [x_async,x_async_comparedtoOptimalOverTime] = parallelFunctionForSims_computations(labindex,c,A_assignment,Aineq_local,beq,bineq,EIPS,regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
end
% delete(gcp);

x_async_1 = cell2mat(x_async(1,1));
x_async_0_75 = cell2mat(x_async(1,2));
x_async_0_5 = cell2mat(x_async(1,3));
x_async_0_1 = cell2mat(x_async(1,4));

xComparedToOptimalOverTime_1 = cell2mat(x_async_comparedtoOptimalOverTime(1,1));
x_async_comparedtoOptimalOverTime_0_75 = cell2mat(x_async_comparedtoOptimalOverTime(1,2));
x_async_comparedtoOptimalOverTime_0_5 = cell2mat(x_async_comparedtoOptimalOverTime(1,3));
x_async_comparedtoOptimalOverTime_0_1 = cell2mat(x_async_comparedtoOptimalOverTime(1,4));


save('parallelHundredvaryingCompRates')

% [x_centralized, lamGradTest,xComparedToOptimalOverTime,lamList] = centralized_PD_Alg(xInitial,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);
if all(Aineq_local*x_async_0_1<=bineq) && all(A_assignment*x_async_0_1<=beq)
    disp("Async Form is Feasible comp 0.1"); 
end
if all(Aineq_local*x_async_0_75<=bineq) && all(A_assignment*x_async_0_75<=beq)
    disp("Async Form is Feasible comp 0.75");
end
if all(Aineq_local*x_async_0_5<=bineq)&& all(A_assignment*x_async_0_5<=beq)
    disp("Async Form is Feasible comp 0.5");
end
if all(Aineq_local*x_async_1<=bineq) && all(A_assignment*x_async_1<=beq)
    disp("Centralized Form is Feasible");
end
%%
%% Parallel Tests spmd tests
%Varying slaters values
scBlock = false;
primalNum = robotNum;
dualNum = robotNum-30;
commRate = 1.0;
plotVals = false;
compute_rate = 1.0;
% p = parpool(4);
spmd(4)
%     labindex
% % spmdIndex = 1;
   [x_async,x_async_comparedtoOptimalOverTime] = parallelFunctionForSims_slaters(labindex,c,A_assignment,Aineq_local,beq,bineq,EIPS,regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
end
% delete(gcp);

x_async_1 = cell2mat(x_async(1,1));
x_async_0_75 = cell2mat(x_async(1,2));
x_async_0_5 = cell2mat(x_async(1,3));
x_async_0_1 = cell2mat(x_async(1,4));

xComparedToOptimalOverTime_1 = cell2mat(x_async_comparedtoOptimalOverTime(1,1));
x_async_comparedtoOptimalOverTime_0_75 = cell2mat(x_async_comparedtoOptimalOverTime(1,2));
x_async_comparedtoOptimalOverTime_0_5 = cell2mat(x_async_comparedtoOptimalOverTime(1,3));
x_async_comparedtoOptimalOverTime_0_1 = cell2mat(x_async_comparedtoOptimalOverTime(1,4));


save('parallelHundredvaryingSlaterVals')

% [x_centralized, lamGradTest,xComparedToOptimalOverTime,lamList] = centralized_PD_Alg(xInitial,lamInit,c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,xIntegerOptimal,primalReg,dualReg);
if all(Aineq_local*x_async_0_1<=bineq) && all(A_assignment*x_async_0_1<=beq)
    disp("Slater's term 3 is Feasible"); 
end
if all(Aineq_local*x_async_0_75<=bineq) && all(A_assignment*x_async_0_75<=beq)
    disp("Slater's term 1 is Feasible");
end
if all(Aineq_local*x_async_0_5<=bineq)&& all(A_assignment*x_async_0_5<=beq)
    disp("Slater's term 2 is Feasible");
end
if all(Aineq_local*x_async_1<=bineq) && all(A_assignment*x_async_1<=beq)
    disp("Slater's term 0 is Feasible");
end

% %%
% 
% diff_centralized = norm(x_centralized-xIntegerOptimal,1)/(1e-10+norm(xIntegerOptimal,1));
% diff_async2_central = norm(x_centralized-x_async);
% 
% 
% 
% %%
% %Suboptimality Error Bound
% % Hoffman(Aconst)
% h = 1;
% % kap = sqrt(p*q);
% kap = 1;
% apriori_optimalityEB = (norm(c,2)*(Hoffman([A_assignment;Aineq_local])*norm([A_assignment;Aineq_local],inf)+kap)+norm(c,2)*kap)*(h/2);
% apriori_feasibilityEB = (Hoffman([A_assignment;Aineq_local])*norm([A_assignment;Aineq_local],inf)+2*kap)*(h/2);

%A posteriori assumes we know the optimal solution centralized
% %Note: I am taking the 2-norm for distance instead of the infimum
% apost_optimalityEB_central = norm(c,2)*(norm(x_centralized-xIntegerOptimal,2)+kap*0.5)+norm(c,2)*kap*(h/2);
% apost_feasibilityEB_central = norm(x_centralized-xIntegerOptimal,2)+kap*h;

%Async form
% apost_optimalityEB_async = norm(c,2)*(norm(x_async_0_1-xIntegerOptimal,2)+kap*0.5)+norm(c,2)*kap*(h/2);
% apost_feasibilityEB_async = norm(x_async_0_1-xIntegerOptimal,2)+kap*h;
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



% %%
% close all
% figure ()
% fontSizeOverall = 25;
% hold on
% plotLen = 5000;
% if length(xComparedToOptimalOverTime_1)>= plotLen
%     semilogy(xComparedToOptimalOverTime_1(2:plotLen)','-','LineWidth',2)
%     semilogy(x_async_comparedtoOptimalOverTime_0_75(2:plotLen)','-','LineWidth',2)
%     semilogy(x_async_comparedtoOptimalOverTime_0_5(2:plotLen)','-','LineWidth',2)
%     semilogy(x_async_comparedtoOptimalOverTime_0_1(2:plotLen)','-','LineWidth',2)
% 
% else
%     semilogy(xComparedToOptimalOverTime_1(2:end)','-','LineWidth',2)
%     semilogy(x_async_comparedtoOptimalOverTime_0_75(2:end)','-','LineWidth',2)
%     semilogy(x_async_comparedtoOptimalOverTime_0_5(2:end)','-','LineWidth',2)
%     semilogy(x_async_comparedtoOptimalOverTime_0_1(2:end)','-','LineWidth',2)
% end
% 
% %Plot Theoretical Error Bounds
% % EB_plot_length = length(xComparedToOptimalOverTime_1(2:end));
% % eb_xPlot = linspace(0,EB_plot_length,EB_plot_length);
% % eb_yPlot = apost_feasibilityEB_async *ones(EB_plot_length);
% % plot(eb_xPlot,eb_yPlot,'g')
% 
% title('Distance to Optimum','FontSize',fontSizeOverall)
% ylabel('$||z_{\kappa}(k) - z^{*}_{MILP}||$','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
% xlabel('Iteration Number','FontWeight','Bold','FontSize',fontSizeOverall)
% ax = gca;
% ax.FontSize = fontSizeOverall;
% % ylabel('Distance Between Iterates','FontWeight','Bold')
% % hold off
% legend('Comm Rate=1.0','Comm Rate=0.75','Comm Rate=0.5','Comm Rate=0.1','A-Post Feasibility','Interpreter','latex','FontSize',fontSizeOverall);
% % legend('Comm Rate=1.0','Comm Rate=0.1','Interpreter','latex','FontSize',fontSizeOverall);

% %%
% %Testing varying computation rate
% fontSizeOverall = 25;
% % close all
% %Async Primal Dual Algorithm Parameters
% scBlock = false;
% primalNum = robotNum;
% dualNum = robotNum;
% commRate = 0.5;
% plotVals = false;
% compute_rate = 0.75;
% maxIterations = 10^5;
% tolerance = 10^-5;
% %Async Primal Dual Algorithm
% [x_async, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal);
% 
% all(Aineq_local*x_async<=bineq) && all(A_assignment*x_async<=beq)
% figure()
% semilogy(x_async_comparedtoOptimalOverTime(2:end)','-','LineWidth',2)
% titleFull = strcat('Distance to Optimum with Comm Rate: ',num2str(commRate),', and Comp Rate: ',num2str(compute_rate));
% title(titleFull,'FontSize',fontSizeOverall)
% ylabel('$||z_{\kappa}(k) - z^{*}_{MILP}||$','Interpreter','latex','FontWeight','Bold','FontSize',fontSizeOverall)
% xlabel('Iteration Number','FontWeight','Bold','FontSize',fontSizeOverall)