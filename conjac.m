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
    x1 = X(ix1);
    u1 = X(iu1);
    
    J(ic,ix1) = eye(4);
    J(ic(end)+1,iu1) = 1;
    ic = ic+nstates+1;
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        ix2 = ix1 + nvarpernode;
		iu2 = iu1 + nvarpernode;
		x1 = X(ix1);
		x2 = X(ix2);
		u1 = X(iu1);
		u2 = X(iu2);
        
        omega = params.omega(:,NperSU*(j-1)+i);
        [f, dfdx, dfdxdot, dfdu] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);
        
        J(ic,ix1) = dfdx/2 - dfdxdot/h;
        J(ic,ix2) = dfdx/2 + dfdxdot/h;
        J(ic,iu1) = dfdu/2;
        J(ic,iu2) = dfdu/2;
        ix1 = ix1+nvarpernode;
        iu1 = iu1+nvarpernode;
        ic = ic+nstates;
    end
    
    %Last node should have pi/2 and zero and ??
    x1 = X(ix1);
    u1 = X(iu1);
    J(ic,ix1) = [1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0];
    J(ic(end)+1,iu1) = 1;
    
    ix1 = ix1+nvarpernode;
    iu1 = iu1+nvarpernode;
    ic = ic+nstates;
end
