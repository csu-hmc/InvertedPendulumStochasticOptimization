function [c] = confun(X, params)
% define constraints

NperSU = params.NperSU;
NSU = params.NSU;
nstates = params.nstates;
optstates = params.optstates;
ncontrols = params.ncontrols;
nvarperSU = params.nvarperSU;
h = params.h;

ix = 1:optstates;
iu = optstates+(1:ncontrols);
ic = 1:nstates;
        
nvarpernode = params.nvarpernode;
ncon = params.ncon;

c = zeros(ncon,1);
% Constraints on dynamics, only BE so far
for j = 1:NSU
    % First node should match initial conditions
    x1 = X(ix);
    u1 = X(iu);
    c(ic) = [x1(1)+pi/2; x1(2)];
    c(ic(end)+(1:ncontrols)) = u1;
    
    if j == 1
        ic = ic+nstates+ncontrols;
    else
        x2 = X(ix-nvarperSU);
        c(ic(end)+ncontrols+(1:(optstates-nstates))) = [x1(3)-x2(3);x1(4)-x2(4)];
        ic = ic+optstates+ncontrols;
    end    
    
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        x1 = X(ix);
		u1 = X(iu);
		x2 = X(ix+nvarpernode);
		u2 = X(iu+nvarpernode);
      
        omega = params.omega(:,NperSU*(j-1)+i);
        dyns = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);
        c(ic) = dyns;
        
        %Constraints for desired states, first SU follows dynamics, all others
        %match first
        if j > 1
            u2 = X(iu-nvarperSU);
            x2 = X(ix-nvarperSU);
            c(ic(end)+ncontrols+(1:(optstates-nstates))) = [x1(3)-x2(3);x1(4)-x2(4)];
            c(ic(end)+(1:ncontrols)) = u1 - u2;
            ic = ic+optstates-nstates+ncontrols;
        end      
        ix = ix+nvarpernode;
        iu = iu+nvarpernode;
        ic = ic+nstates;
    end
    
    %Average of last nodes should equal pi/2, u should be zero
    x1 = X(ix);
    u1 = X(iu);
    c(end) = c(end)+1/NSU*x1(1);
    c(ic(1):ic(ncontrols)) = u1;
    
    if j > 1
        x2 = X(ix-nvarperSU);        
        c(ic(ncontrols)+(1:(optstates-nstates))) = [x1(3)-x2(3);x1(4)-x2(4)];
        ic = ic+optstates-nstates;
    end
 
    ix = ix+nvarpernode;
    iu = iu+nvarpernode;
    ic = ic+ncontrols;         
end

c(end) = c(end)-pi/2;

