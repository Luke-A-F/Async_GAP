function [x,fval,exitflag,info] = opti_scip(H,f,A,rl,ru,lb,ub,xint,sos,qc,x0,opts)
%OPTI_SCIP Solve a LP/MILP/QP/MIQP/QCQP/MIQCQP using SCIP
%
%   min 0.5*x'*H*x + f'*x      subject to:     rl <= A*x <= ru
%    x                                         qrl <= x'Qx + l'x <= qru
%                                              lb <= x <= ub
%                                              for i = 1..n: xi in Z
%                                              for j = 1..m: xj in {0,1}
%
%   x = opti_scip([],f,A,rl,ru,lb,ub,xint,x0) solves a LP/MILP where f is the
%   objective vector, A,rl,ru are the linear constraints, lb,ub are the
%   bounds and xint is a string of integer variables ('C', 'I', 'B'). The
%   point x0 provides a primal solution.
%
%   x = opti_scip(H,f,A,rl,ru,lb,ub,xint,x0) solves a QP/MIQP where H is the
%   objective matrix, and the remainder of the arguments are as above.
%
%   x = opti_scip(H,...,xint,sos,x0) sos is a structure with fields type,
%   index, and weight for SOS.
%
%   x = opti_scip(H,...,sos,qc,x0) qc is structure with fields Q, l, and qrl
%   and qru for quadratic constraints.
%
%   x = opti_scip(H,f,...,qc,x0,opts) uses opts to pass optiset options to the
%   solver.
%
%   [x,fval,exitflag,info] = opti_scip(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   This is a wrapper for SCIP using the mex interface.
%   See the license for SCIP.
%   Copyright (C) 2012 Jonathan Currie (IPL)

t = tic;

% handle missing arguments
if nargin < 12, opts = optiset; end
if nargin < 11, x0 = []; end
if nargin < 10, qc = []; end
if nargin < 9, sos = []; end
if nargin < 8, xint = repmat('C',size(f)); end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, error('You must supply at least 5 arguments to opti_scip.'); end

warn = strcmpi(opts.warnings,'all');

% check sparsity
if(~issparse(A))
    if(warn)
        optiwarn('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)].');
    end
    A = sparse(A);
end
if(~isempty(H) && ~issparse(H))
    if(warn)
        optiwarn('opti:sparse','The H matrix should be sparse, correcting: [sparse(A)].');
    end
    H = sparse(H);
end
if(~isempty(qc))
    if(~isfield(qc,'Q') || ~isfield(qc,'l') || ~isfield(qc,'qrl') || ~isfield(qc,'qru'))
        error('The qc structure must contain fields Q, l, qrl and qru.');
    end
    if(~isempty(qc.Q))
        err = 0;
        if(iscell(qc.Q))
            for i=1:length(qc.Q)
                if(~issparse(qc.Q{i}))
                    err = 1;
                    qc.Q{i} = sparse(qc.Q{i});
                end
            end
        elseif(~issparse(qc.Q))
            err = 1;
            qc.Q = sparse(qc.Q);
        end
        if(err && warn)
            optiwarn('opti:sparse','The Q matrix should be sparse, correcting: [sparse(Q)].');
        end
    end
end

% remove H if all nz
if(~isempty(H) && nnz(H) == 0), H = []; end

% MEX contains error checking
[x,fval,exitflag,stats] = scip(H, f, A, rl, ru, lb, ub, xint, sos, qc, [], x0, opts);

% assign outputs
info.BBNodes = stats.BBnodes;
info.BBGap = stats.BBgap;
info.PrimalBound = stats.PrimalBound;
info.DualBound = stats.DualBound;
info.Time = toc(t);
info.Algorithm = 'SCIP: Spatial Branch and Bound';

% process return code
[info.Status,exitflag] = scipRetCode(exitflag);

% return display level
function print_level = dispLevel(lev)
switch(lower(lev))
    case 'off'
        print_level = 0;
    case 'iter'
        print_level = 4;
    case 'final'
        print_level = 3;
end
