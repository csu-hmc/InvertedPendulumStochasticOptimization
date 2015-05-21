function result = runOpt(X0, L, U, params)

options = optimoptions(@fmincon,'Algorithm','sqp','display', 'iter');

[x,fval,exitflag] = fmincon(@(X) objfun(X,params), X0, [], [], [], [], L, U, @(X) confun(X,params), options);
result.X = x;
result.fval = fval;
result.exitflag = exitflag;

% funcs.objective = @(X) objfun(X,params);
% funcs.gradient  = @(X) objgrad(X,params);
% funcs.constraints = @(X) confun(X, params);
% funcs.jacobian    = @(X) conjac(X, params);
% funcs.jacobianstructure = @() conjacstructure(L,U,params);
% options.lb = L;
% options.ub = U;
% options.cl = zeros(params.ncon,1);
% options.cu = zeros(params.ncon,1);	
% options.ipopt.max_iter = 10000;
% options.ipopt.hessian_approximation = 'limited-memory';
% options.ipopt.mu_strategy = 'adaptive';		% worked better than 'monotone'
% % options.ipopt.tol = problem.Tol;
% options.ipopt.linear_solver = 'mumps';
% % options.ipopt.constr_viol_tol = problem.ConstraintTol;
% % options.ipopt.compl_inf_tol = problem.ConstraintTol;
% options.ipopt.print_level = 0;
% options.ipopt.bound_frac = 0.001;			% worked better than 0.01 or 0.0001
% options.ipopt.bound_push = options.ipopt.bound_frac;
% [X, info] = ipopt(X0,funcs,options);
% disp(['IPOPT status: ' (info.status)]);
% result.info = info.status;
% result.X = X;