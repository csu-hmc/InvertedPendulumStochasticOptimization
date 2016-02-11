function result = Optimize(X0, L, U, params)

if strcmp(params.optimizer, 'IPOPT')
    funcs.objective = @(X) objfun(X,params);
    funcs.gradient  = @(X) objgrad(X,params);
    funcs.constraints = @(X) confun(X, params);
    funcs.jacobian    = @(X) conjac(X, params);
    funcs.jacobianstructure = @() conjacstructure(L,U,params);
    options.lb = L;
    options.ub = U;
    [cl,cu] = getconstraintsb(params);
    options.cl = cl;%zeros(params.ncon,1);
    options.cu = cu;%zeros(params.ncon,1);	
    options.ipopt.max_iter = 5000;
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy = 'adaptive';		% worked better than 'monotone'
    options.ipopt.tol = 1e-4;
    options.ipopt.linear_solver = 'mumps'; %'ma57';% 
    options.ipopt.constr_viol_tol = 1e-3;
    options.ipopt.compl_inf_tol = 1e-3;
    options.ipopt.print_level = 5;
    options.ipopt.bound_frac = 0.01;%1e-8;
    options.ipopt.bound_push = options.ipopt.bound_frac;
    options.ipopt.recalc_y = 'yes';
    options.ipopt.dual_inf_tol = 1e-1;
    options.ipopt.compl_inf_tol = 1e-2;
    % options.ipopt.bound_relax_factor = 0;
    [X, info] = ipopt(X0,funcs,options);
    disp(['IPOPT status: ' num2str(info.status)]);
    result.info = info.status;
    result.X = X;
    result.params = params;
    result.obj = objfun(X, params);
elseif strcmp(params.optimizer, 'SNOPT')
    testspec.spc = which('testspec.spc');
    snspec ( testspec.spc );
    % Output informative files
    snprint   ([params.snoptname '.out']);
    snsummary ([params.snoptname '.sum']);
    [cl,cu] = getconstraintsb(params);
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
elseif strcmp(params.optimizer, 'WORHP')
    result = theworhpthing(X0,L,U,params);
else
    error('Idiot, I do not recognize this solver')
end