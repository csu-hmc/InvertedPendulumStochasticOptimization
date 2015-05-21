function result = Optimize(X0, L, U, params)

funcs.objective = @(X) objfun(X,params);
funcs.gradient  = @(X) objgrad(X,params);
funcs.constraints = @(X) confun(X, params);
funcs.jacobian    = @(X) conjac(X, params);
funcs.jacobianstructure = @() conjacstructure(L,U,params);
options.lb = L;
options.ub = U;
options.cl = zeros(params.ncon,1);
options.cu = zeros(params.ncon,1);	
for i = 1:params.NSU
    if i == 1
        options.cl(params.nconSU1-params.nvarpernode+1) = -0.1;
        options.cu(params.nconSU1-params.nvarpernode+1) = 0.1;
%         options.cl(params.nconSU1-params.nvarpernode+2) = -0.3;
%         options.cu(params.nconSU1-params.nvarpernode+2) = 0.3;
    else
        options.cl(params.nconSU1+(i-1)*params.nconperSU-params.nvarpernode+1) = -0.1;
        options.cu(params.nconSU1+(i-1)*params.nconperSU-params.nvarpernode+1) = 0.1;
%         options.cl(params.nconSU1+(i-1)*params.nconperSU-params.nvarpernode+2) = -0.3;
%         options.cu(params.nconSU1+(i-1)*params.nconperSU-params.nvarpernode+2) = 0.3;
    end
end
options.ipopt.max_iter = 500000;
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
options.ipopt.dual_inf_tol = 1e-2;
options.ipopt.compl_inf_tol = 1e-2;
% options.ipopt.bound_relax_factor = 0;
[X, info] = ipopt(X0,funcs,options);
disp(['IPOPT status: ' num2str(info.status)]);
result.info = info.status;
result.X = X;
result.params = params;