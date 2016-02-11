function [c, Bu,Bu2] = confun(X, params)
% define constraints

NperSU = params.NperSU;
NSU = params.NSU;
nstates = params.nstates;
ncontrols = params.ncontrols;    
nvarpernode1= params.nvarpernode1;
nvarpernode = params.nvarpernode;
nconSU = params.nconSU;
nvarSU1 = params.nvarSU1;
nvarSU = params.nvarSU;
h = params.h;
omega = params.omega;

c = zeros(params.ncon,1);
% Constraints on dynamics, only BE so far
for j = 1:NSU
    if j == 1
        ix = (1:nstates);
    else
        ix = nvarSU1+(j-2)*nvarSU+(1:nstates);
    end
    % First node should match initial conditions
    x1 = X(ix);
    ic = nconSU*(j-1)+(1:nstates);
    c(ic) = [x1(1)+pi/2; x1(2)];
    
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        if j == 1
            ix = (i-1)*nvarpernode1+(1:nstates);
        else
            ix = nvarSU1+(j-2)*nvarSU+(i-1)*nvarpernode+(1:nstates);
        end
        iu = (i-1)*nvarpernode1+nstates+(1:ncontrols);
        iKs = (i-1)*nvarpernode1+nstates+ncontrols+(1:2);
        ic = nconSU*(j-1)+i*nstates+(1:nstates);
        
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
      
        omega_now = omega(:,NperSU*(j-1)+i);
        u1 = findTorque(u01, K1, x1);
        u2 = findTorque(u02, K2, x2);
        [dyns,Bu(:,i)] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega_now, params);
        c(ic) = dyns;
    end
    
    if j == 1
        ix = (NperSU-1)*nvarpernode1+(1:nstates);
    else
        ix = nvarSU1+(j-2)*nvarSU+(NperSU-1)*nvarpernode+(1:nstates);
    end
    
    %Average of last nodes should equal pi/2
    x1 = X(ix);
    c(end-1) = c(end-1)+1/NSU*x1(1);
    % Velocity should be zero on average
    c(end) = c(end)+1/NSU*x1(end);
end

c(end-1) = c(end-1)-pi/2;
params.obj1 = params.obj;
% norm(c)

