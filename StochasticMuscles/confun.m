function [c] = confun(X, params)
% define constraints

NperSU = params.NperSU;
NSU = params.NSU;
nstates = params.nstates;
ncontrols = params.ncontrols;
h = params.h;

ix = 1:nstates;
iu = nstates+(1:ncontrols);
ic = 1:nstates;
icineq = params.nconeq+(1:params.nmus);
        
nvarpernode = params.nvarpernode;
nvarpernode1 = params.nvarpernode1;
ncon = params.ncon;
ndof = params.ndof;

c = zeros(ncon,1);

% Constraints on dynamics, only Midpoint so far
for j = 1:NSU
    % First node should match initial conditions
    x1 = X(ix);
    c(ic(1:ndof*2)) = [x1(1)-pi/2; x1(2)];%3*
    ic = ic+ndof*2; 
    
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        x1 = X(ix);
		u1 = X(iu);
        if j == 1
            x2 = X(ix+nvarpernode1);
        else
            x2 = X(ix+nvarpernode);
        end
		u2 = X(iu+nvarpernode1);
        
        u1 = findTorque(u1,X(end-1),X(end),x1(1:ndof*2),params);
        u2 = findTorque(u2,X(end-1),X(end),x2(1:ndof*2),params);
      
        omega = params.omega(:,NperSU*(j-1)+i);
        dyns = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);
        c(ic) = dyns;
        c(icineq) = u1;
        
        if j == 1
            ix = ix+nvarpernode1;
        else
            ix = ix+nvarpernode;
        end
        iu = iu+nvarpernode1;
        ic = ic+nstates;
        icineq = icineq+params.nmus;
    end
    
    %Average of last nodes should equal pi/2
    x1 = X(ix);
    u1 = X(iu);
    u1 = findTorque(u1,X(end-1),X(end),x1(1:ndof*2),params);
    
    c(icineq) = u1;
%     c(params.nconeq-1) = c(params.nconeq-1)+1/NSU*x1(1);
%     % Velocity should be zero on average
%     c(params.nconeq) = c(params.nconeq)+1/NSU*x1(2);
%  
    if j == 1
        ix = ix+nvarpernode1;
    else
        ix = ix+nvarpernode;
    end
    iu = nstates+(1:ncontrols);
    icineq = icineq+params.nmus;
end

% c(params.nconeq-1) = c(params.nconeq-1)-pi/2;

