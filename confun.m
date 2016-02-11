function [c] = confun(X, params)
% define constraints

NperSU = params.NperSU;
NSU = params.NSU;
nstates = params.nstates;
ncontrols = params.ncontrols;
nvarperSU = params.nvarperSU;
h = params.h;

ix = 1:nstates;
iu = nstates+(1:ncontrols);
ic = 1:nstates;
        
nvarpernode = params.nvarpernode;
ncon = params.ncon;

c = zeros(ncon,1);
% Constraints on dynamics, only Midpoint so far
for j = 1:NSU
    % First node should match initial conditions
    x1 = X(ix);
    c(ic) = [x1(1)+pi/2; x1(2)];%3*
    ic = ic+nstates; 
    
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        x1 = X(ix);
		u1 = X(iu);
		x2 = X(ix+nvarpernode);
		u2 = X(iu+nvarpernode);
      
        omega = params.omega(:,NperSU*(j-1)+i);
        dyns = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params, X(end-1), X(end));
        c(ic) = dyns;
              
        %Open loop controls should match
        if j > 1
            u2 = X(iu-nvarperSU);
            c(ic(end)+(1:ncontrols)) = u1-u2;
            ic = ic+ncontrols;
        end
        
        ix = ix+nvarpernode;
        iu = iu+nvarpernode;
        ic = ic+nstates;
    end
    
    %Average of last nodes should equal pi/2
    x1 = X(ix);
    u1 = X(iu);
    c(end-1) = c(end-1)+1/NSU*x1(1);
    % Velocity should be zero on average
    c(end) = c(end)+1/NSU*x1(end);
    if j > 1
        % Control should be equal to previous
        u2 = X(iu-nvarperSU);
        c(ic(1:ncontrols)) = u1-u2;
        ic = ic+ncontrols;
    end
    
    %inequality constraint
    if params.ineq == 1
        c(ic(1):ic(1)+params.NperSU-1) = ineqconfun(X, params,j);
        ic = ic+params.NperSU;
    end
 
    ix = ix+nvarpernode;
    iu = iu+nvarpernode;
end

c(end-1) = c(end-1)-pi/2;

