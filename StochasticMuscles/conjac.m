function J = conjac(X, params)
% define constraints jacobian

NperSU = params.NperSU;
NSU = params.NSU;
h = params.h; %time step

nvarpernode = params.nvarpernode;
nvarpernode1 = params.nvarpernode1;
nvars = params.nvars;
ncon = params.ncon;
nstates = params.nstates;
ncontrols = params.ncontrols;
ndof = params.ndof;

ix1 = 1:nstates;
iu1 = nstates+(1:ncontrols);
ic = 1:nstates;
icineq = params.nconeq+(1:params.nmus);

J = spalloc(ncon,nvars,params.Jnnz);

% Constraints on dynamics, only BE so far
for j = 1:NSU
    
    %First node should be at initial condition
    J(ic(1:ndof*2),ix1(1:ndof*2)) = eye(ndof*2);
    ic = ic+ndof*2;
    
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        iu2 = iu1 + nvarpernode1;
        if j == 1
            ix2 = ix1 + nvarpernode1;
        else
            ix2 = ix1 + nvarpernode;
        end
        x1 = X(ix1);
        x2 = X(ix2);
        u1 = X(iu1);
        u2 = X(iu2);

        [u1, du1dx, du1dK, du1dKd] =  findTorque(u1,X(end-1),X(end),x1(1:ndof*2),params);
        [u2, du2dx, du2dK, du2dKd] =  findTorque(u2,X(end-1),X(end),x2(1:ndof*2),params);

        omega = params.omega(:,NperSU*(j-1)+i);
        [f, dfdx, dfdxdot, dfdu] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);

        dfdK = dfdu/2*du1dK+dfdu/2*du2dK;
        dfdKd= dfdu/2*du1dKd+dfdu/2*du2dKd;
                
        J(ic,ix1) = dfdx/2 - dfdxdot/h + dfdu/2*[du1dx zeros(params.nmus,params.nmus*2)];
        J(ic,ix2) = dfdx/2 + dfdxdot/h + dfdu/2*[du2dx zeros(params.nmus,params.nmus*2)];
        J(ic,iu1) = dfdu/2;
        J(ic,iu2) = dfdu/2;
        J(ic,end-1) = dfdK;
        J(ic,end) = dfdKd;
        
        %inequality constraint on u1
        J(icineq,iu1) = eye(ncontrols);
        J(icineq,ix1(1:ndof*2)) = du1dx;
        J(icineq,end-1) = du1dK;
        J(icineq,end) = du1dKd;

        ix1 = ix2;
        iu1 = iu2;
        ic = ic+nstates;
        icineq = icineq+params.nmus;
    end
    
    x1 = X(ix1);
    u1 = X(iu1);
    [~, du1dx, du1dK, du1dKd] =  findTorque(u1,X(end-1),X(end),x1(1:ndof*2),params);
    %inequality constraint on u1
    J(icineq,iu1) = eye(ncontrols);
    J(icineq,ix1(1:ndof*2)) = du1dx;
    J(icineq,end-1) = du1dK;
    J(icineq,end) = du1dKd;
        
    %Average of last nodes should equal pi/2, u should be zero
%     J(params.nconeq-1,ix1(1)) = 1/NSU;
%     J(params.nconeq,ix1(2)) = 1/NSU;
%     
    if j == 1
        ix1 = ix1+nvarpernode1;
    else
        ix1 = ix1+nvarpernode;
    end
    iu1 = nstates+(1:ncontrols);
    icineq = icineq+params.nmus;
end


