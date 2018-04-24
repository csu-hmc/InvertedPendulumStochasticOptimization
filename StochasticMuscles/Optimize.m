function result = Optimize(X0, L, U, params)

if strcmp(params.solver, 'IPOPT')
    funcs.objective = @(X) objfun(X,params);
    funcs.gradient  = @(X) objgrad(X,params);
    funcs.constraints = @(X) confun(X, params);
    funcs.jacobian    = @(X) conjac(X, params);
    funcs.jacobianstructure = @() conjacstructure(params);
    options.lb = L;
    options.ub = U;
    options.cl = zeros(params.ncon,1);
    options.cu = zeros(params.ncon,1);
%     options.cu(end) = 0.1;
%     options.cl(end) = -0.1;
    options.ipopt.max_iter = 1000;%00;
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy = 'adaptive';
    options.ipopt.tol = 1e-2;
    options.ipopt.linear_solver = 'mumps'; 
    options.ipopt.constr_viol_tol = 1e-4;
    options.ipopt.dual_inf_tol = 1e-2;
    options.ipopt.compl_inf_tol = 1e-4;
    options.ipopt.print_level = 5;
    options.ipopt.bound_frac = 0.01;
    options.ipopt.bound_push = options.ipopt.bound_frac;
    if params.warmstart
        options.ipopt.warm_start_init_point = 'yes';
        options.ipopt.warm_start_bound_frac = 1e-16;
        options.ipopt.warm_start_bound_push = 1e-16;
        options.ipopt.warm_start_mult_bound_push = 1e-16;
        options.ipopt.warm_start_slack_bound_frac = 1e-16;
        options.ipopt.warm_start_slack_bound_push = 1e-16;
        options.zl = params.zl;
        options.zu = params.zu;
        options.lambda = params.lambda;
    end
    [X, info] = ipopt(X0,funcs,options);
    disp(['IPOPT status: ' num2str(info.status)]);
    result.info = info.status;
    result.zl = info.zl;
    result.zu = info.zu;
    result.lambda = info.lambda;
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
    cu = zeros(params.ncon,1);
    cl(end) = -.1;
    cu(end) = .1;

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