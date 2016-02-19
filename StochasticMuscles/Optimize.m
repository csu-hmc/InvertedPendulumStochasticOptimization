function result = Optimize(X0, L, U, params)

if strcmp(params.solver, 'IPOPT')
    funcs.objective = @(X) objfun(X,params);
    funcs.gradient  = @(X) objgrad(X,params);
    funcs.constraints = @(X) confun(X, params);
    funcs.jacobian    = @(X) conjac(X, params);
    funcs.jacobianstructure = @() conjacstructure(L,U,params);
    options.lb = L;
    options.ub = U;
    options.cl = zeros(params.ncon,1);
    options.cu = [zeros(params.nconeq,1); 2*ones(params.nconineq,1)];	
    options.ipopt.max_iter = 10000;
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy = 'adaptive';
    options.ipopt.tol = 1e-4;
    options.ipopt.linear_solver = 'mumps'; 
    options.ipopt.constr_viol_tol = 1e-3;
    options.ipopt.compl_inf_tol = 1e-3;
    options.ipopt.print_level = 5;
    options.ipopt.bound_frac = 0.01;
    options.ipopt.bound_push = options.ipopt.bound_frac;
    options.ipopt.recalc_y = 'yes';
    options.ipopt.dual_inf_tol = 1e-1;
    options.ipopt.compl_inf_tol = 1e-2;
    [X, info] = ipopt(X0,funcs,options);
    disp(['IPOPT status: ' num2str(info.status)]);
    result.info = info.status;
    result.X = X;
    result.params = params;
    result.obj = objfun(X, params);
else
    Prob = conAssign(@(X) objfun(X,params), @(X) objgrad(X,params), [], [], L, U, 'pendulum', X0, [], 0, ...
                [], [], [], @(X) confun(X, params), @(X) conjac(X, params), [], params.Jpattern, ...
                zeros(params.ncon,1), zeros(params.ncon,1), ...
                [], [], [],[]);
    % Prob.SOL.optPar(1)= 1;		% uncomment this to get snoptsum.txt and snoptpri.txt
    Prob.SOL.optPar(9) = 1e-3;		% feasibility tolerance
    Prob.SOL.optPar(10) = 1e-4;		% optimality tolerance
    Prob.SOL.optPar(11) = 1e-3;     % Minor feasibility tolerance (1e-6)
    Prob.SOL.optPar(12) = 1e-4;		% optimality tolerance
    Prob.SOL.optPar(30) = 1000000;  % maximal sum of minor iterations (max(10000,20*m))
    Prob.SOL.optPar(35) = 2000;
    Prob.SOL.optPar(36) = 40000; % maximal number of minor iterations in the solution of the QP problem (500)
    Prob.SOL.moremem = 10000000; % increase internal memory
    Result = tomRun('snopt',Prob,3);
    X = Result.x_k;
    disp(Result.ExitText);
    result.info = Result.ExitFlag;
    result.X = X;
    result.params = params;
    result.obj = objfun(X, params);
end

%     if params.ineq == 1
%         for i = 1:params.NSU
%             if i == 1
%                 options.cl((i-1)*(params.nstates*params.NperSU+params.ncontrols+params.NperSU)+params.nstates*params.NperSU+params.ncontrols-1+(1:params.NperSU)) =  -0.1;%inf;
%                 options.cu((i-1)*(params.nstates*params.NperSU+params.ncontrols+params.NperSU)+params.nstates*params.NperSU+params.ncontrols-1+(1:params.NperSU)) = 0.1;%inf;
%             else
%                 options.cl((i-1)*((params.nstates+params.ncontrols)*params.NperSU+params.NperSU)+params.nstates*params.NperSU+params.ncontrols-1+(1:params.NperSU)) = -0.1;%inf;
%                 options.cu((i-1)*((params.nstates+params.ncontrols)*params.NperSU+params.NperSU)+params.nstates*params.NperSU+params.ncontrols-1+(1:params.NperSU)) = 0.1;%inf;
%             end
%         end
%     end