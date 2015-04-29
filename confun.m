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
        
nvarpernode = params.nvarpernode;
ncon = params.ncon;

c = zeros(ncon,1);
% Constraints on dynamics, only BE so far
for j = 1:NSU
    % First node should match initial conditions
    x1 = X(ix);
    u1 = X(iu);
    
    c(ic) = [x1(1); x1(2); x1(3); x1(4)-10];
    c(ic(end)+1) = u1;
    ic = ic+nstates+1;
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        x1 = X(ix);
		u1 = X(iu);
		x2 = X(ix+nvarpernode);
		u2 = X(iu+nvarpernode);
      
        omega = params.omega(:,NperSU*(j-1)+i);
        c(ic) = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);
        
        ix = ix+nvarpernode;
        iu = iu+nvarpernode;
        ic = ic+nstates;
    end
    
    %Last node should have pi/2 and zero and ??
    x1 = X(ix);
    u1 = X(iu);
    c(ic) = [x1(1)-pi/2;x1(2);0;0];
    c(ic(end)+1) = u1;
    
    ix = ix+nvarpernode;
    iu = iu+nvarpernode;
    ic = ic+nstates;
end

