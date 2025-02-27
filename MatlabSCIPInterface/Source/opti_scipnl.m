function [x,fval,exitflag,info] = opti_scipnl(fun,A,rl,ru,lb,ub,nlcon,cl,cu,xint,x0,xval,opts)
%OPTI_SCIPNL Solve a NLP/MINLP using SCIP to Global Optimality
%
%   min fun(x)                 subject to:     rl <= A*x <= ru
%    x                                         lb <= x <= ub
%                                              cl <= nlcon(x) <= cu
%                                              for i = 1..n: xi in Z
%                                              for j = 1..m: xj in {0,1} [i!=j]
%
%   Full Calling Form:
%     [x,fval,exitflag,info] = opti_scipnl(fun,A,rl,ru,lb,ub,nlcon,cl,cu,xint,x0,xval,opts)
%
%   NOTE:
%     This solver (and interface) parses your supplied function(s) into an
%     algebraic description of the problem. This means your functions must be
%     deterministic (no random numbers or stochastic elements), and must also
%     contain only a subset of functions (log,log10,exp,abs,sqrt,norm). In
%     addition, only MATLAB code may be parsed, meaning no Simulink or MEX
%     models may be optimized using this solver.
%
%
%   x = opti_scipnl(fun) solves an unconstrained NLP (UNO) with objective
%   specified by fun. Note you must supply xval if solving an UNO.
%
%   x = opti_scipnl(fun,A,rl,ru) solves a NLP subject to linear constraints
%   where A is the linear constraint matrix, and rl and ru are lower and
%   upper row bounds.
%
%   x = opti_scipnl(fun,A,rl,ru,lb,ub) solves subject to upper and lower
%   decision variable bounds. Finite bounds on ALL decision variables are
%   HIGHLY recommended when solving global optimization problems.
%
%   x = opti_scipnl(fun,...,ub,nlcon,cl,cu) solves subject to nonlinear
%   constraints specified by nlcon (column vector of constraint functions)
%   and lower and upper row bounds (cl, cu).
%
%   x = opti_scipnl(fun,...,cu,xint) solves a MINLP subject to integer and
%   binary constraints specified by xint. Integer constraints are specified
%   using a string of integer variables ('C', 'I', 'B').
%
%   x = opti_scipnl(fun,...,xint,x0) uses x0 as a starting point.
%
%   x = opti_scipnl(fun,...,xint,xval) uses xval as a validation
%   vector of the decision variables for the MEX interface, to use to
%   verify your objective and constraints have been successfully
%   converted to SCIP expressions.
%
%   x = opti_scipnl(fun,...,x0,opts) uses opts to pass optiset options to the
%   solver.
%
%   [x,fval,exitflag,info] = opti_scipnl(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   This is a wrapper for SCIP using the mex interface.
%   See the SCIP license.
%   Copyright (C) 2012/2013 Jonathan Currie (IPL)

t = tic;

% Handle missing arguments
if nargin < 13, opts = optiset; end
if nargin < 12, xval = []; end
if nargin < 11, x0 = []; end
if nargin < 10, xint = []; end
if nargin < 9, cu = []; end
if nargin < 8, cl = []; end
if nargin < 7, nlcon = []; end
if nargin < 6, ub = []; end
if nargin < 5, lb = []; end
if nargin < 4, ru = []; end
if nargin < 3, rl = []; end
if nargin < 2, A = []; end
if nargin < 1, error('You must supply at least one argument to opti_scipnl.'); end

% check for number of decision variables
if(~isempty(xval))       % vector or matrix
    ndec = numel(xval);
elseif(~isempty(lb))   % vector
    ndec = length(lb);
elseif(~isempty(ub))   % vector
    ndec = length(ub);
elseif(~isempty(xint)) % vector
    ndec = length(xint);
elseif(~isempty(A))
    ndec = size(A,2);
else
    error('The MATLAB-SCIP Interface requires the number of variables to be specified.\n%s',...
          ['The interface cannot determine this number from your problem thus you will need',...
          ' to specify this via xval, e.g. xval = ones(no_vars,1).']);
end

% build xval if empty, used for validation
if(isempty(xval))
    xval = randn(ndec,1);
end

if(isfield(opts,'warnings'))
    warn = strcmpi(opts.warnings,'all');
else
    warn = 0;
end

% check for NaNs
if(any(isnan(xval)))
    if(warn)
        optiwarn('opti:nan','For the SCIP interface xval must not contain NaN. Replacing with random numbers.');
    end
    xval = randn(ndec,1);
end

% encode nonlinear constraints into SCIP MEX interface instruction lists
x = scipvar(size(xval));
if(~isempty(nlcon))
    try
        n = nlcon(x);
    catch ME
        ex = processEqErr(ME,'a constraint');
        throwAsCaller(ex);
    end
    if(numel(n) > 1) % multiple constraints
        nl.instr = cell(numel(n),1);
        for i = 1:numel(n)
            if(isnumeric(n(i))) % read as number
                nl.instr{i} = [0 n(i)];
            elseif(isempty(n(i).ins)) % assume just a variable
                nl.instr{i} = [1; n(i).indx];
            else
                nl.instr{i} = n(i).ins;
            end
        end
    else
        if(isnumeric(n)) % read as number
            nl.instr = [0 n];
        elseif(isempty(n.ins)) % assume just a variable
            nl.instr = [1; n.indx];
        else
            nl.instr = n.ins;
        end
    end
    nl.cl = cl;
    nl.cu = cu;
    % verification fields
    nl.nlcon_val = nlcon(xval);
    nl.xval = xval;
    % check for Inf or NaN
    if(any(isnan(nl.nlcon_val)) || any(isinf(nl.nlcon_val)))
        error('One or more constraints resulted in Inf or NaN at the initial guess (xval). Please provide a better initial guess vector.');
    end
end

% encode nonlinear objective into SCIP MEX interface instruction list
if(~isempty(fun))
    try
        f = fun(x);
    catch ME
        ex = processEqErr(ME,'the objective');
        throwAsCaller(ex);
    end
    if(isnumeric(f)) %read as number
        nl.obj_instr = [0 f];
    elseif(isempty(f.ins)) %assume just a variable
        nl.obj_instr = [1; f.indx];
    else
        nl.obj_instr = f.ins;
    end
    % verification field
    nl.obj_val = fun(xval);
    nl.xval = xval;
    % check for Inf or NaN
    if(any(isnan(nl.obj_val)) || any(isinf(nl.obj_val)))
        error('The objective resulted in Inf or NaN at the initial guess (xval). Please provide a better initial guess vector.');
    end
end

% check sparsity
if(~isempty(A) && ~issparse(A))
    if(warn)
        optiwarn('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)].');
    end
    A = sparse(A);
end

% adding SCIP settings if specified
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    sopts = scipset(opts.solverOpts);
else
    sopts = [];
end

% add OPTI options
if(isfield(opts,'maxtime') && ~isempty(opts.maxtime))
    sopts.maxtime = opts.maxtime;
end
if(isfield(opts,'maxiter') && ~isempty(opts.maxiter))
    sopts.maxiter = opts.maxiter;
end
if(isfield(opts,'maxnodes') && ~isempty(opts.maxnodes))
    sopts.maxnodes = opts.maxnodes;
end
if(isfield(opts,'tolrfun') && ~isempty(opts.tolrfun))
    sopts.tolrfun = opts.tolrfun;
end
if(isfield(opts,'objbias') && ~isempty(opts.objbias))
    sopts.objbias = opts.objbias;
end
if(isfield(opts,'display') && ~isempty(opts.display))
    sopts.display = dispLevel(opts.display);
end
sopts.optiver = optiver;

% run SCIP
[x,fval,exitflag,stats] = scip([],zeros(ndec,1),A,rl,ru,lb,ub,xint,[],[],nl,x0,sopts);

% reshape output
x = reshape(x,size(xval));

% assign outputs
info.BBNodes = stats.BBnodes;
info.BBGap = stats.BBgap;
info.PrimalBound = stats.PrimalBound;
info.DualBound = stats.DualBound;
info.Time = toc(t);
info.Algorithm = 'SCIP: Spatial Branch and Bound';

% process return code
[info.Status,exitflag] = scipRetCode(exitflag);


% return CLP compatible display level
function  print_level = dispLevel(lev)
switch(lower(lev))
    case'off'
        print_level = 0;
    case 'iter'
        print_level = 4;
    case 'final'
        print_level = 3;
end

function ex = processEqErr(ME,name)

str = [];
for i = 1:length(ME.stack)
    if(i==1)
        if(length(ME.stack)>1 && strcmpi(ME.stack(2).name,'opti_scipnl'))
            str = sprintf('%sError using "%s" (line %d):\n%s\n',str,ME.stack(i).name,ME.stack(i).line,ME.message);
        else
            str = sprintf('%sError using "%s" (line %d):\n%s\n\nError Trace:\n',str,ME.stack(i).name,ME.stack(i).line,ME.message);
        end
    else
        if(strcmpi(ME.stack(i).name,'opti_scipnl'))
            break;
        end
        str = sprintf('%s- Error in "%s" (line %d)\n',str,ME.stack(i).name,ME.stack(i).line);
    end
end

ex = MException('OPTI:SCIPVAR',['There was an error processing %s function into SCIP compatible form.\n'...
                                'Please examine the below error to correct your function:\n\n%s\nRemember '...
                                'SCIP only supports a subset of MATLAB commands. Use:\n>> methods(scipvar)\n\nto see compatible functions and operators.\n\n'],name,str);
