function result = Optimize(X0, L, U, params)

if strcmp(params.solver, 'IPOPT')
    funcs.objective = @(X) objfun(X,params);
    funcs.gradient  = @(X) objgrad(X,params);
    funcs.constraints = @(X) confun(X, params);
    funcs.jacobian    = @(X) conjac(X, params);
    funcs.jacobianstructure = @() conjacstructure(params);
%     funcs.hessian = @(X,sigma,lambda) objhess(X,sigma,lambda,params);
%     funcs.hessianstructure = @() hessstructure(L, U, params);
    options.lb = L;
    options.ub = U;
    options.cl = zeros(params.ncon,1);
    options.cu = zeros(params.ncon,1);
    options.cl(end) = -.1;
%     if params.asat
%         options.cu = [zeros(params.nconeq,1); 100*ones(params.nconineq,1)];	
%     else
%         options.cu = [zeros(params.nconeq,1); 1*ones(params.nconineq,1)];	
%     end
    options.cu(end) = .1;
    options.ipopt.max_iter = 2000;%00;
    options.ipopt.hessian_approximation = 'limited-memory';
%     options.ipopt.limited_memory_max_history = 100;
%     options.ipopt.limited_memory_max_skipping = 1;
%     options.ipopt.recalc_y = 'no';
%     options.ipopt.first_hessian_perturbation = 1e-9;
    options.ipopt.mu_strategy = 'adaptive';
    options.ipopt.tol = 1e-2;
    options.ipopt.linear_solver = 'mumps'; 
    options.ipopt.constr_viol_tol = 1e-5;
    options.ipopt.dual_inf_tol = 1e-2;
    options.ipopt.compl_inf_tol = 1e-5;
    options.ipopt.print_level = 5;
    options.ipopt.bound_frac = 0.01;
    options.ipopt.bound_push = options.ipopt.bound_frac;
%     options.ipopt.slack_bound_frac = 1e-9;
%     options.ipopt.slack_bound_push = options.ipopt.slack_bound_frac;
%     options.ipopt.mu_init = 1e-9;
    [X, info] = ipopt(X0,funcs,options);
    disp(['IPOPT status: ' num2str(info.status)]);
    result.info = info.status;
    result.X = X;
    result.params = params;
    result.obj = objfun(X, params);
else
% elseif strcmp(params.solver, 'SNOPT')
    testspec.spc = which('testspec.spc');
    snspec ( testspec.spc );
    % Output informative files
    snprint   ([params.snoptname '.out']);
    snsummary ([params.snoptname '.sum']);
    cl = zeros(params.ncon,1);
    if params.asat
        cu = [zeros(params.nconeq,1); 100*ones(params.nconineq,1)];	
    else
        cu = [zeros(params.nconeq,1); 1*ones(params.nconineq,1)];	
    end
    FL = [-inf;cl];
    FU = [inf;cu];
    xmul = zeros(size(L));
    Fmul = zeros(size(FL));
    xstate = zeros(size(X0));
    Fstate = zeros(size(FL));
    A = [];
    iAfun = [];
    jAvar = [];
    iGfun = params.DGrow;
    jGvar = params.DGcol;
    if params.warmstart == 1;
        snset ('Warm start')
    end
    [X,F,INFO] = snopt(X0,L,U,xmul, xstate,FL,FU,Fmul, Fstate, @objconfunp, ...
        A, iAfun, jAvar, iGfun, jGvar );
    snprint   off;
    snsummary off;
    result.info = INFO;
    result.obj = F(1);
    result.X = X;
    result.params = params;
end