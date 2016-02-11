function J = conjac_eq(X, params)

% define constraints jacobian

NperSU = params.NperSU;
NSU = params.NSU;
nstates = params.nstates;
ncontrols = params.ncontrols;    
nvarpernode1= params.nvarpernode1;
nvarpernode = params.nvarpernode;
h = params.h;

ix1 = 1:nstates;
iu1 = nstates+(1:ncontrols);
iKs1 = nstates+ncontrols+(1:2);
ic = 1:nstates;

J = spalloc(params.ncon,params.nvars,params.Jnnz);

% Constraints on dynamics, only ME so far
for j = 1:NSU
    
    %First node should be at initial condition
    J(ic,ix1) = eye(nstates);
    ic = ic+nstates;
    
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        if j == 1
            ix2 = ix1 + nvarpernode1;         
        else
            ix2 = ix1 + nvarpernode;
        end
        iu2 = iu1 + nvarpernode1;
        iKs2= iKs1+ nvarpernode1;
        x1 = X(ix1);
		u01 = X(iu1);
		x2 = X(ix2);
		u02 = X(iu2);
        K1 = X(iKs1);
        K2 = X(iKs2);

        omega = params.omega(:,NperSU*(j-1)+i);
        [u1, du1dK1, du1dx1] = findTorque(u01, K1, x1);
        [u2, du2dK2, du2dx2] = findTorque(u02, K2, x2);
        [~, ~,dfdx, dfdxdot, dfdu,dfdK] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);

        J(ic,ix1) = dfdx/2 - dfdxdot/h + dfdu/2*du1dx1;
        J(ic,ix2) = dfdx/2 + dfdxdot/h + dfdu/2*du2dx2;
        J(ic,iu1) = dfdu/2;
        J(ic,iu2) = dfdu/2;
        J(ic,iKs1) = dfdK/2 + dfdu/2*du1dK1;
        J(ic,iKs2) = dfdK/2 + dfdu/2*du2dK2;
        
        if j == 1
            ix1 = ix1+nvarpernode1;
        else
            ix1 = ix1+nvarpernode;
        end
        iu1 = iu1+nvarpernode1;
        iKs1= iKs1+nvarpernode1;
        ic = ic+nstates;
    end
       
    %Last node should equal pi/2 and 0 velocity
    J(ic,ix1) = eye(nstates);
    ic = ic + nstates;
    
    if j == 1
        ix1 = ix1+nvarpernode1;
    else
        ix1 = ix1+nvarpernode;
    end
    %Reset controls and Ks index
    iu1 = nstates+(1:ncontrols);
    iKs1 = nstates+ncontrols+(1:2);
end
% 
% 
