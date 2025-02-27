%% OPTI Toolbox Linear Programming Examples
%
% This file contains a number of LP problems and demonstrates how to solve
% them using the OPTI Toolbox.
%
%   Copyright (C) 2014 Jonathan Currie (IPL)

%% Example 1
% This is a simple two decision variable LP which will use for the next few
% examples.
clc
% problem
f = -[6 5]';                % objective function (min f'x)
A = [1,4; 6,4; 2, -5];      % linear inequality constraints (Ax <= b)
b = [16;28;6];
lb = [0;0];                 % bounds on x (lb <= x <= ub)
ub = [10;10];

% Problems are built using the opti constructor:
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub)

% Call solve to solve the problem:
[x,fval,exitflag,info] = solve(Opt)


%% Example 2 - Mixed Constraints
% There are multiple ways to enter linear constraints. This includes by
% groups such as inequalities and equalities (1), as row linear constraints
% (2), mixed constraints (3) or individually (4). All variants will produce
% an equivalent model for the optimizer.
clc
% problem
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
Aeq = [1 1 1];
beq = 40;
lb = [0;0;0];
ub = [40;inf;inf];

% 1)
Opt = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub)

% 2)
Ar = [A;Aeq];
rl = [-Inf(size(b));beq];
ru = [b;beq];
Opt = opti('f',f,'lin',Ar,rl,ru,'bounds',lb,ub)

% 3)
Amix = [A;Aeq];
bmix = [b;beq];
e = [-1;-1;0];  % specify via constraint type vector (-1 <=, 0 ==, 1 >=)
Opt = opti('f',f,'mix',Amix,bmix,e,'bounds',lb,ub)

% 4)
Opt = opti('f',f,'A',A,'b',b,'Aeq',Aeq,'beq',beq,'bounds',lb,ub)

%% Example 3 - Maximization
% You can use the 'sense' option to change between minimization and
% maximization. Use sense = -1 for maximization, 1 for minimization.
clc
fmax = -f; % identical but with inverted sign for this example

Opt = opti('f',fmax,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'sense',-1)

x = solve(Opt)

%% Example 4 - Sparse LPs
% All solvers are setup to directly solve sparse systems, which is the
% preferred format for most solvers:
clc
% A larger sparse LP for the next few examples:
load sparseLP1;

opts = optiset('solver','scip'); % Solve with SCIP
Opt = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'options',opts);
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% Example 5 - Using the Optimization Toolbox
% If you have the MATLAB optimization toolbox you can also solve OPTI
% problems using it:
clc
opts = optiset('solver','scip');
Opt = opti(Opt,opts)
[x,fval,exitflag,info] = solve(Opt);
fval
info
