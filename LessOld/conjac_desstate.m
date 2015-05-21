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
    J(ic(1):ic(end)+ncontrols,ix1(1):iu1(end)) = eye(nstates+ncontrols);
    ic = ic+nstates+ncontrols;
    
    % Dynamics should match next node till one before last node
    for i = 1:NperSU-1
        ix2 = ix1 + nvarpernode;
        iu2 = iu1 + nvarpernode;
        x1 = X(ix1);
        x2 = X(ix2);
        u1 = X(iu1);
        u2 = X(iu2);

        omega = params.omega(:,NperSU*(j-1)+i);
        [f, dfdx, dfdxdot, dfdu, dfdK,dfdKd] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, X(end-1), X(end), omega, params);

        J(ic(1:nstates/2),ix1) = dfdx(1:nstates/2,:)/2 - dfdxdot(1:nstates/2,:)/h;
        J(ic(1:nstates/2),ix2) = dfdx(1:nstates/2,:)/2 + dfdxdot(1:nstates/2,:)/h;
        J(ic(1:nstates/2),iu1) = dfdu(1:nstates/2,:)/2;
        J(ic(1:nstates/2),iu2) = dfdu(1:nstates/2,:)/2;
        J(ic(1:nstates/2),end-1) = dfdK(1:nstates/2);
        J(ic(1:nstates/2),end) = dfdKd(1:nstates/2);
        
        %Constraints for desired states, first SU follows dynamics, all others
        %match first
        if j == 1                    
            J(ic(nstates/2+1:end),ix1) = dfdx(nstates/2+1:end,:)/2 - dfdxdot(nstates/2+1:end,:)/h;
            J(ic(nstates/2+1:end),ix2) = dfdx(nstates/2+1:end,:)/2 + dfdxdot(nstates/2+1:end,:)/h;
            J(ic(nstates/2+1:end),iu1) = dfdu(nstates/2+1:end,:)/2;
            J(ic(nstates/2+1:end),iu2) = dfdu(nstates/2+1:end,:)/2;
            J(ic(nstates/2+1:end),end-1) = dfdK(nstates/2+1:end);
            J(ic(nstates/2+1:end),end) = dfdKd(nstates/2+1:end);
            ic = ic+nstates;
        else
            ix = ix1(nstates/2+1:end);
            iu = iu1(1:ncontrols/2);
            
            J(ic(nstates/2+1:end), ix) = eye(nstates/2);
            J(ic(nstates/2+1:end), ix-nvarperSU) = -eye(nstates/2);
            J(ic(end)+ncontrols/2,iu) = 1;
            J(ic(end)+ncontrols/2,iu-nvarperSU) = -1;
            ic = ic+nstates+ncontrols/2;
        end      
        ix1 = ix1+nvarpernode;
        iu1 = iu1+nvarpernode;
    end
    
    %Last node should have pi/2 and zero and ?? and zero input?
    if j == 1
        J(ic(1:nstates),ix1) = [eye(2), zeros(2);zeros(2) zeros(2)]; %zeros(4);%
        J(ic(nstates)+(1:ncontrols/2),iu1(1:ncontrols/2)) = eye(ncontrols/2);
    else
        J(ic(1:nstates), ix1) = eye(4); %[zeros(2) zeros(2); zeros(2) eye(2)];%
        J(ic(nstates/2+1:end),ix1-nvarperSU) = [zeros(nstates/2) -eye(nstates/2)];
        J(ic(nstates)+(1:ncontrols),iu1) = [eye(ncontrols/2) zeros(ncontrols/2); eye(ncontrols/2) zeros(ncontrols/2)];
        J(ic(nstates)+(ncontrols/2+1:ncontrols),iu1-nvarperSU)= [-eye(ncontrols/2) zeros(ncontrols/2)];
    end
    
    ix1 = ix1+nvarpernode;
    iu1 = iu1+nvarpernode;
    ic = ic+nstates+ncontrols;
end
