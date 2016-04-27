function J = conjac(X, params)
% keyboard
% define constraints jacobian

h = params.h; %time step
N = params.N;

nvarpernode = params.nvarpernode;
nvars = params.nvars;
ncon = params.ncon;
nstates = params.nstates;
ncontrols = params.ncontrols;
ndof = params.ndof;

ix1 = ncontrols+(1:nstates);
iu = 1:ncontrols;
u = X(iu);
ic = 1:nstates;
icineq = params.nconeq+(1:params.nmus);

J = spalloc(ncon,nvars,params.Jnnz);

for j = 1:params.Nsamples
    % Constraints on dynamics
    % First node should be at initial condition
    J(ic,ix1) = eye(nstates); %(1:ndof*2)
%     J(ic(3:4),iu) = -ones(params.nmus,ncontrols);
    ic = ic+nstates;

    % Dynamics should match next node till one before last node
    for i = 1:N-1
        if i == N
            ix2 = ncontrols+(1:nstates);
        else
            ix2 = ix1 + nvarpernode;
        end

        x1 = X(ix1);
        x2 = X(ix2);

        [u1, du1dx, du1du, du1dK, du1dKd] =  findTorque(u,X(end-params.Ks+1:end),x1(1:ndof*2),params);
        [u2, du2dx, du2du, du2dK, du2dKd] =  findTorque(u,X(end-params.Ks+1:end),x2(1:ndof*2),params);

        omega = params.omega(:,(j-1)*params.N+i);

        if strcmp(params.method, 'ME')
            [~, dfdx, dfdxdot, dfdu] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);
            dfdu1 = dfdu/2;
            dfdu2 = dfdu/2;
            dfdK = dfdu1/2*du1dK+dfdu2/2*du2dK;
            dfdKd= dfdu1/2*du1dKd+dfdu2/2*du2dKd;


            J(ic,ix1) = dfdx/2 - dfdxdot/h + dfdu/2*[du1dx zeros(params.nmus,params.nmus*2)];
            J(ic,ix2) = dfdx/2 + dfdxdot/h + dfdu/2*[du2dx zeros(params.nmus,params.nmus*2)];
        else
            [~, dfdx, dfdxdot, dfdu] = StocDyn(x2,(x2-x1)/h,u2, omega, params);
            dfdu1 = zeros(size(dfdu));
            dfdu2 = dfdu;
            dfdK = dfdu2*du2dK;
            dfdKd= dfdu2*du2dKd;
            
            J(ic,ix1) = -dfdxdot/h;
            J(ic,ix2) = dfdx + dfdxdot/h + dfdu*[du2dx zeros(params.nmus,params.nmus*2)];
        end

        if params.ncontrols == 1
            J(ic,iu) = (dfdu1*du1du+dfdu2*du2du)*[1;1];
        else
            J(ic,iu) = (dfdu1*du1du+dfdu2*du2du);
        end
        J(ic,end-params.Ks+1:end) = [dfdK dfdKd];

%         J(end,ix1(1)) = 1/N/params.Nsamples;
%         %inequality constraint on u1
%         J(icineq,iu) = eye(ncontrols);
%         J(icineq,ix1(1:ndof*2)) = du1dx;
%         J(icineq,end-1) = du1dK;
%         J(icineq,end) = du1dKd;

        ix1 = ix2;
        ic = ic+nstates;
%         icineq = icineq+params.nmus;
    end

    x1 = X(ix1);
%     [~, du1dx, du1dK, du1dKd] =  findTorque(u,X(end-1:end),x1(1:ndof*2),params);
    %inequality constraint on u1
%     J(icineq,iu) = eye(ncontrols);
%     J(icineq,ix1(1:ndof*2)) = du1dx;
%     J(icineq,end-1) = du1dK;
%     J(icineq,end) = du1dKd;
%     J(end,ix1(1)) = 1/params.Nsamples/N; %
%     icineq = icineq+params.nmus;
    ix1 = ix1+nvarpernode;
end
% 
J(end,iu) = [u(2)/u(1)^2 -1/u(1)];