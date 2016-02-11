function [c, Bu,Bu2] = confun_eq(X, params)
% define constraints

NperSU = params.NperSU;
NSU = params.NSU;
nstates = params.nstates;
ncontrols = params.ncontrols;    
nvarpernode1= params.nvarpernode1;
nvarpernode = params.nvarpernode;
h = params.h;

ix = 1:nstates;
iu = nstates+(1:ncontrols);
iKs = nstates+ncontrols+(1:2);
ic = 1:nstates;

c = zeros(params.ncon,1);
% Constraints on dynamics, only BE so far
for j = 1:NSU
    % First node should match initial conditions
    x1 = X(ix);
    c(ic) = [x1(1)+pi/2; x1(2)];
    ic = ic+nstates; 
    
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        x1 = X(ix);
		u01 = X(iu);
        K1 = X(iKs);
        u02 = X(iu+nvarpernode1);
        K2 = X(iKs+nvarpernode1);
        if j == 1
            x2 = X(ix+nvarpernode1);
        else
            x2 = X(ix+nvarpernode);
        end
      
        omega = params.omega(:,NperSU*(j-1)+i);
        u1 = findTorque(u01, K1, x1);
        u2 = findTorque(u02, K2, x2);
        [dyns,Bu(:,i)] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);
        c(ic) = dyns;
        
        if j == 1
            ix = ix+nvarpernode1;
        else
            ix = ix+nvarpernode;
        end
        iu = iu+nvarpernode1;
        iKs= iKs+nvarpernode1;
        ic = ic+nstates;
    end
    
    %Last node should equal pi/2 and 0 velocity
    x1 = X(ix);
    c(ic) = x1-[pi/2;0];
    ic = ic+nstates;

    if j == 1
        ix = ix+nvarpernode1;
    else
        ix = ix+nvarpernode;
    end
    %Reset controls and Ks index
    iu = nstates+(1:ncontrols);
    iKs = nstates+ncontrols+(1:2);
end

