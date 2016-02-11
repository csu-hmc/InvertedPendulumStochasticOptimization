function J = conjac(X, params)

% define constraints jacobian
NperSU = params.NperSU;
NSU = params.NSU;
nstates = params.nstates;
ncontrols = params.ncontrols;    
nvarpernode1= params.nvarpernode1;
nvarpernode = params.nvarpernode;
nconSU = params.nconSU;
nconeq = params.nconeq;
nvarSU1 = params.nvarSU1;
nvarSU = params.nvarSU;
h = params.h;
omega = params.omega;

J = spalloc(params.ncon,params.nvars,params.Jnnz);

% Constraints on dynamics, only ME so far
for j = 1:NSU
    %First node should be at initial condition
    if j == 1
        ix1 = (1:nstates);
    else
        ix1 = nvarSU1+(j-2)*nvarSU+(1:nstates);
    end
    ic = nconSU*(j-1)+(1:nstates);
    J(ic,ix1) = eye(nstates);
    
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        if j == 1
            ix1 = (i-1)*nvarpernode1+(1:nstates);
            ix2 = ix1 + nvarpernode1;
        else
            ix1 = nvarSU1+(j-2)*nvarSU+(i-1)*nvarpernode+(1:nstates);
            ix2 = ix1 + nvarpernode;
        end
        iu1 = (i-1)*nvarpernode1+nstates+(1:ncontrols);
        iKs1 = (i-1)*nvarpernode1+nstates+ncontrols+(1:2);
        ic = nconSU*(j-1)+i*nstates+(1:nstates);
        iu2 = iu1 + nvarpernode1;
        iKs2= iKs1+ nvarpernode1;
        x1 = X(ix1);
		u01 = X(iu1);
		x2 = X(ix2);
		u02 = X(iu2);
        K1 = X(iKs1);
        K2 = X(iKs2);

        omega_now = omega(:,NperSU*(j-1)+i);
        [u1, du1dK1, du1dx1] = findTorque(u01, K1, x1);
        [u2, du2dK2, du2dx2] = findTorque(u02, K2, x2);
        [~, ~,dfdx, dfdxdot, dfdu,dfdK] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega_now, params);

        J(ic, ix1) = dfdx/2-dfdxdot/h+dfdu/2*du1dx1;
        J(ic, ix2) = dfdx/2+dfdxdot/h+dfdu/2*du2dx2;
        J(ic, iu1) = dfdu/2;
        J(ic, iu2) = dfdu/2;
        J(ic, iKs1)= dfdK/2+dfdu/2*du1dK1;
        J(ic, iKs2)= dfdK/2+dfdu/2*du2dK2;
        
    end
    
    if j == 1
        ix1 = (NperSU-1)*nvarpernode1+(1:nstates);
    else
        ix1 = nvarSU1+(j-2)*nvarSU+(NperSU-1)*nvarpernode+(1:nstates);
    end
    
    %Average of last nodes should equal pi/2, u should be zero
    J(nconeq-1,ix1(1)) = 1/NSU;
    J(nconeq,ix1(2)) = 1/NSU;
    
    %Inequality constraint on the last node
    J(nconeq+j,ix1(1)) = 1;
end
% 
% 
