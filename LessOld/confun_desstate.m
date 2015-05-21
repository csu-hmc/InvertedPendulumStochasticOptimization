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
% Constraints on dynamics, only BE so far
for j = 1:NSU
    % First node should match initial conditions
    x1 = X(ix);
    u1 = X(iu);
    
    c(ic) = [x1(1); x1(2); x1(3); x1(4)];
    c(ic(end)+(1:ncontrols)) = u1;
    ic = ic+nstates+ncontrols;
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        x1 = X(ix(1:nstates));
		u1 = X(iu(1:ncontrols));
		x2 = X(ix(1:nstates)+nvarpernode);
		u2 = X(iu(1:ncontrols)+nvarpernode);
      
        omega = params.omega(:,NperSU*(j-1)+i);
        dyns = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, X(end-1), X(end), omega, params);
        c(ic(1:nstates/2)) = dyns(1:nstates/2);
        
        %Constraints for desired states, first SU follows dynamics, all others
        %match first
        if j == 1
            c(ic(nstates/2+1:end)) = dyns(nstates/2+1:end);
            ic = ic+nstates;
        else
            x1 = X(ix(nstates/2+1:end));
            u1 = X(iu(1:ncontrols/2));
            x2 = X(ix(nstates/2+1:end)-nvarperSU);
            u2 = X(iu(1:ncontrols/2)-nvarperSU);
            c(ic(nstates/2+1:end)) = x1 - x2;
            c(ic(end)+1) = u1 - u2;
            ic = ic+nstates+ncontrols/2;
        end      
        ix = ix+nvarpernode;
        iu = iu+nvarpernode;
    end
    
    %Last node should have pi/2 and zero and ?? and zero input?
    if j == 1
        x1 = X(ix(1:nstates));
        u1 = X(iu(1:ncontrols));
        c(ic(1:nstates)) = [x1(1)-pi/2;x1(2);0;0];%x1(3)-pi/2;x1(4)];% [0;0;0;0];% 
        c(ic(nstates)+(1:ncontrols/2)) = u1(1:ncontrols/2);
    else
        x1 = X(ix(1:nstates/2));
        u1 = X(iu(1:ncontrols/2));
        x1des = X(ix(nstates/2+1:end));
        x2des = X(ix(nstates/2+1:end)-nvarperSU);
        u1des = X(iu(1:ncontrols/2));
        u2des = X(iu(1:ncontrols/2)-nvarperSU);
        c(ic(1:nstates)) = [x1(1)-pi/2;x1(2);x1des-x2des]; % [0;0;x1des-x2des];% 
        c(ic(nstates)+(1:ncontrols)) = [u1(1:ncontrols/2);u1des-u2des];
    end
        
    ix = ix+nvarpernode;
    iu = iu+nvarpernode;
    ic = ic+(nstates+ncontrols);         
end

