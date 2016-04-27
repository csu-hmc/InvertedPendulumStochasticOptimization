function J = conjac_tdfeed(X, params)
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
    % Constraints on dynamics, only mE so far
    % First node should be at initial condition
    J(ic,ix1) = eye(nstates); %(1:ndof*2)
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

        if i == 1
            ixp = ix1;
            [u1, du1dx, du1dK, du1dKd] = findTorque(u,X(end-1:end),x1(1:ndof*2),params);
            [u2, du2dx, du2dK, du2dKd] = findTorque(u,X(end-1:end),x1(1:ndof*2),params);
        else
            ixp = ix1 - nvarpernode;
            xprev = X(ixp);
            [u1, du1dx, du1dK, du1dKd] = findTorque(u,X(end-1:end),xprev(1:ndof*2),params);
            [u2, du2dx, du2dK, du2dKd] = findTorque(u,X(end-1:end),x1(1:ndof*2),params);
        end

        omega = params.omega(:,(j-1)*params.N+i);

        if strcmp(params.method, 'ME')
            [~, dfdx, dfdxdot, dfdu] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);
            dfdK = dfdu/2*du1dK+dfdu/2*du2dK;
            dfdKd= dfdu/2*du1dKd+dfdu/2*du2dKd;

            if i == 1
                J(ic,ix1) = dfdx/2 - dfdxdot/h + dfdu*[du2dx+du1dx zeros(params.nmus,params.nmus*2)];
                J(ic,ix2) = dfdx/2 + dfdxdot/h;
            else
                J(ic,ix1) = dfdx/2 - dfdxdot/h + dfdu/2*[du2dx zeros(params.nmus,params.nmus*2)];
                J(ic,ix2) = dfdx/2 + dfdxdot/h;
                J(ic,ixp) = dfdu/2*[du1dx zeros(params.nmus,params.nmus*2)];
            end
        else
            [~, dfdx, dfdxdot, dfdu] = StocDyn(x2,(x2-x1)/h,u2, omega, params);
            dfdK = dfdu*du2dK;
            dfdKd= dfdu*du2dKd;

            J(ic,ix1) = -dfdxdot/h + dfdu*[du2dx zeros(params.nmus,params.nmus*2)];
            J(ic,ix2) = dfdx + dfdxdot/h;
        end

        J(ic,iu) = dfdu*[1;1];
        J(ic,end-1) = dfdK;
        J(ic,end) = dfdKd;

        J(end,ix1(1)) = 1/N/params.Nsamples;
        %inequality constraint on u1
        J(icineq,iu) = eye(ncontrols);
        J(icineq,ix1(1:ndof*2)) = du2dx;
        J(icineq,end-1) = du2dK;
        J(icineq,end) = du2dKd;

        ix1 = ix2;
        ic = ic+nstates;
        icineq = icineq+params.nmus;
    end

    x1 = X(ix1);
    [~, du1dx, du1dK, du1dKd] =  findTorque(u,X(end-1:end),x1(1:ndof*2),params);
    %inequality constraint on u1
    J(icineq,iu) = eye(ncontrols);
    J(icineq,ix1(1:ndof*2)) = du1dx;
    J(icineq,end-1) = du1dK;
    J(icineq,end) = du1dKd;
    J(end,ix1(1)) = 1/params.Nsamples/N; %
    icineq = icineq+params.nmus;
    ix1 = ix1+nvarpernode;
end
% 
