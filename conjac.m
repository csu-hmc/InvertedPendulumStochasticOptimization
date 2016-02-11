function J = conjac(X, params)
% define constraints jacobian

NperSU = params.NperSU;
NSU = params.NSU;
h = params.h; %time step

nvarpernode = params.nvarpernode;
nvars = params.nvars;
ncon = params.ncon;
nstates = params.nstates;
ncontrols = params.ncontrols;
nvarperSU = params.nvarperSU;

ix1 = 1:nstates;
iu1 = nstates+(1:ncontrols);
ic = 1:nstates;

J = spalloc(ncon,nvars,params.Jnnz);

% Constraints on dynamics, only BE so far
for j = 1:NSU
    
    %First node should be at initial condition
    J(ic,ix1) = eye(nstates);
    ic = ic+nstates;
    
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        ix2 = ix1 + nvarpernode;
        iu2 = iu1 + nvarpernode;
        x1 = X(ix1);
        x2 = X(ix2);
        u1 = X(iu1);
        u2 = X(iu2);

        omega = params.omega(:,NperSU*(j-1)+i);
        [f, dfdx, dfdxdot, dfdu, dfdK, dfdKd] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params,X(end-1),X(end));

        J(ic,ix1) = dfdx/2 - dfdxdot/h;
        J(ic,ix2) = dfdx/2 + dfdxdot/h;
        J(ic,iu1) = dfdu/2;
        J(ic,iu2) = dfdu/2;
        J(ic,end-1) = dfdK;
        J(ic,end) = dfdKd;

        %Open loop controls should match
        if j > 1
            iu2 = iu1-nvarperSU;
            J(ic(end)+(1:ncontrols),iu1) = eye(ncontrols);
            J(ic(end)+(1:ncontrols),iu2) = -eye(ncontrols);
            ic = ic+ncontrols;
        end

        ix1 = ix1+nvarpernode;
        iu1 = iu1+nvarpernode;
        ic = ic+nstates;
    end
    
    %Average of last nodes should equal pi/2, u should be zero
    J(end-1,ix1(1)) = 1/NSU;
    J(end,ix1(2)) = 1/NSU;
    if j > 1
	    % Control should be equal to previous
        J(ic(1:ncontrols),iu1) = eye(ncontrols);
        J(ic(1:ncontrols),iu1-nvarperSU) = -eye(ncontrols);
        ic = ic+ncontrols;
    end
    
    %Inequality constraint
    if params.ineq == 1
        J(ic(1):ic(1)+params.NperSU-1,:) = ineqconjac(X, params,j);
        ic = ic+params.NperSU;
    end
    
    ix1 = ix1+nvarpernode;
    iu1 = iu1+nvarpernode;
end


