yalmip('clear')
clear all
close all
clc






curTime = 0;
num_cars = 9;

z = binvar(1,num_cars-1);
t = sdpvar(1,num_cars);
theta = sdpvar(1,num_cars);%Slack variables
s = sdpvar(1,num_cars);%Slack variables


% M = 2000;
M = 200;
tgap1 = 4;%seconds
tgap2 = 7.5;%seconds
t_des = 35+(45-35)*abs(randn(num_cars,1));
w1 = 0.5;
w2 = 0.5;
rng(7)




%Constraint construction

const_cont = [];
const_int = [];

objective = 0;

t0 =curTime;


%Note in full formulation this updates with time
tmin = t0 +0.1+(0.5-0.1)*abs(randn(num_cars-1,1));
tmin = [tmin;0];

for i = 1:num_cars
  

 % objective = objective + w2*theta(i)+w1*s(i);
 %No J2 objective term
objective = objective + w1*s(i);
 
 %Constraints with only continuous variables
 % const_cont = [const_cont, (t(i)-t_des(i))<=theta(i),-(t(i)-t_des(i))<=theta(i), (t(i)-t0)<=s(i), tmin<=t(i)];
 %No J2 obj terms
  const_cont = [const_cont,(t(i)-t0)<=s(i), tmin<=t(i)];
 %No tmin constraint
 % const_cont = [const_cont, (t(i)-t_des(i))<=theta(i),-(t(i)-t_des(i))<=theta(i), (t(i)-t0)<=s(i)];

end
for j = 1:num_cars-1
    if j < num_cars
        k = j+1;
    end
        const_cont=[const_cont,tgap1<=t(j)-t(k)];
         %Constraints with binary variables
        const_int = [const_int,tgap2<=t(j)-t(k)+M*z(j),tgap2<=t(j)-t(k)+M*(1-z(j))];
    
end
constraints = [const_int,const_cont];


ops = sdpsettings('solver','gurobi','usex0',1,'savesolveroutput',1);

%Solves internally with gurobi
solve = optimize(constraints,objective,ops);


%
[exportedModel,recoverymodel] = export(constraints,objective,sdpsettings('solver','intlinprog'));

%Algorithm parameters & initialization
V = exportedModel.A;
W = full(exportedModel.b);
fancy_H = exportedModel.c';
xInitial = zeros(length(exportedModel.ub),1); 

%%
xIntegerOptimal = intlinprog_gurobi(exportedModel.c,exportedModel.intcon,full(exportedModel.A),full(exportedModel.b));
%%
%Algorithm parameters & initialization
Aineq = exportedModel.A;
bineq = full(exportedModel.b);
c_vec = exportedModel.c';

% xi = 1-10^-3;

intConstraintIndex = length(const_int);
binIndex = num_cars-1;

Aineq = full(Aineq);

bineq = full(bineq);
%Relaxation constant close to one
xi = 1-10^-3;

n = size(Aineq,2);
c = c_vec';



%Enlarged Inner Parallel Set Computation
EIPS = -sum(abs(Aineq)')'./2;%+xi;%%%Note this can change dependent on problem formulation


%Algorithm parameters & initialization

lamInit = zeros(size(bineq,1),1);

tolerance = 10^-4;
%Binary Variables
ubPrimal_bin =1+(xi-0.5);
lbPrimal_bin =(0.5-xi);
%Continuous Variables Bounds
ubPrimal_con = 75;
lbPrimal_con = 0;



ubDual = 10;
lbDual = 0;
%Feasible for comm 1.0 with 10^5 iterations
primalReg =10^-3; 
dualReg = 10^-2;






stepSize = [10^-2,10^-4];

lb_granCheck = [];
ub_granCheck = [];
for i = 1:length(c)
    if i< binIndex+1
        lb_granCheck(i) = lbPrimal_bin;
        ub_granCheck(i) = ubPrimal_bin;
    else
        lb_granCheck(i) = lbPrimal_con;
        ub_granCheck(i) = ubPrimal_con;
    end

end

granular_bineq = bineq;
granular_bineq(1:intConstraintIndex) = floor(bineq(1:intConstraintIndex))+EIPS(1:intConstraintIndex);
granular_bineq(intConstraintIndex+1:end) = bineq(intConstraintIndex+1:end)+EIPS(intConstraintIndex+1:end);


isGranularCheck = lb_granCheck';
rowwise_normA = vecnorm(Aineq*isGranularCheck);
[minVal,minrow_index] = min(rowwise_normA);

slaterTerms = (c'*isGranularCheck+(primalReg/2)*norm(isGranularCheck)^(2)*norm(c)*ubPrimal_bin)/(-Aineq(minrow_index,:)*isGranularCheck+bineq(minrow_index)+EIPS(minrow_index));
rowwise_normA = vecnorm(Aineq);
[maxVal,maxrow_index] = max(rowwise_normA);
regCorrectionTest = maxVal*slaterTerms*sqrt(dualReg/(2*primalReg))


%Feasible initial point
% xInitial = [0;0;0;1;1;0;0;1;73;62;54;46;37;25;17;9;1;73;66;55;48;40;29;17;9;1]; 
% granular_bineq = bineq;
% granular_bineq(1:intConstraintIndex) = floor(bineq(1:intConstraintIndex))+EIPS(1:intConstraintIndex)+regCorrectionTest;
% granular_bineq(intConstraintIndex+1:end) = bineq(intConstraintIndex+1:end)+EIPS(intConstraintIndex+1:end)+regCorrectionTest;
% 

%Calculate correct dual upper bound with maths
% isGranularCheck = linprog(zeros(size(c)),Aineq,granular_bineq,[],[],lb_granCheck,ub_granCheck);
% rowwise_normA = vecnorm(Aineq*isGranularCheck);
% [minVal,minrow_index] = min(rowwise_normA);
% 
% slaterTerms = (c'*isGranularCheck+(primalReg/2)*norm(isGranularCheck)^(2)*norm(c)*ubPrimal_bin)/(-Aineq(minrow_index,:)*isGranularCheck+bineq(minrow_index)+EIPS(minrow_index));
% rowwise_normA = vecnorm(Aineq);
% [maxVal,maxrow_index] = max(rowwise_normA);
% dualUB_math = maxVal*slaterTerms*sqrt(dualReg/(2*primalReg))/length(lamInit);

%Async Primal Dual Algorithm Parameters
scBlock = false;


primalNum = 10;
dualNum = 10;

plotVals = false;
%%%Change Compute Rate Here
compute_rate = 1.0;



%Async Primal Dual Algorithm


maxIterations = 5*10^5;
% maxIterations = 2*10^6;
%%%%%Change Comm Rate Here
commRate = 0.5;

[x_async, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime] = Async_PD_traffic(c,Aineq,granular_bineq,stepSize,maxIterations,tolerance,lbPrimal_bin,ubPrimal_bin,lbPrimal_con,ubPrimal_con,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal,binIndex,xInitial,bineq);
x_val = x_async;

if all(Aineq*x_async<=bineq)
    disp("Feasible");
else
    disp("Infeasible");
end

xtilde = x_async;


semilogy(x_async_comparedtoOptimalOverTime)
figure()
semilogy(constr_async)
figure()
semilogy(es_async)
% % % %Assigns intlinprog solution back to YALMIP model
% assign(recover(recoverymodel.used_variables),xtilde);

% x_async_comparedtoOptimalOverTime_01 = x_async_comparedtoOptimalOverTime;
% constr_async_01 = constr_async;
% es_async_01 = es_async;
% save("traffic_MILP_comm_01.mat","x_async_comparedtoOptimalOverTime_01","constr_async_01","es_async_01");
