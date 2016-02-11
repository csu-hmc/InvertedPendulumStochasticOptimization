function J = conjac(X, params)

% define constraints jacobian
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

    if j == 1
        XJ = X(1:nvarSU1);
    else
        XJ = X(nvarSU1+(j-2)*nvarSU+(1:nvarSU));
    end
    omegaJ = omega(:,NperSU*(j-1)+(1:NperSU));
    Xcon = X(1:nvarSU1);
    Jx1 = zeros(nstates,nstates*2, NperSU-1);
    Jx2 = zeros(nstates,nstates*2, NperSU-1);
    Ju1 = zeros(ncontrols,nstates*2+ncontrols, NperSU-1);
    Ju2 = zeros(ncontrols,nstates*2+ncontrols, NperSU-1);
    JK1 = zeros(nstates,nstates*2, NperSU-1);
    JK2 = zeros(nstates,nstates*2, NperSU-1);
    % Dynamics should match next node till one before last node
    parfor i = 1:NperSU-1
        if j == 1
            ix1 = (i-1)*nvarpernode1+(1:nstates);
            ix2 = ix1 + nvarpernode1;
        else
            ix1 = (i-1)*nvarpernode+(1:nstates);
            ix2 = ix1 + nvarpernode;
        end
        iu1 = (i-1)*nvarpernode1+nstates+(1:ncontrols);
        iKs1 = (i-1)*nvarpernode1+nstates+ncontrols+(1:2);
        ic = nconSU*(j-1)+i*nstates+(1:nstates);
        iu2 = iu1 + nvarpernode1;
        iKs2= iKs1+ nvarpernode1;
        x1 = XJ(ix1);
		u01 = Xcon(iu1);
		x2 = XJ(ix2);
		u02 = Xcon(iu2);
        K1 = Xcon(iKs1);
        K2 = Xcon(iKs2);

        omega_now = omegaJ(:,i);
        [u1, du1dK1, du1dx1] = findTorque(u01, K1, x1);
        [u2, du2dK2, du2dx2] = findTorque(u02, K2, x2);
        [~, ~,dfdx, dfdxdot, dfdu,dfdK] = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega_now, params);

        Jx1(:,:,i) = [ic' ix1' dfdx/2-dfdxdot/h+dfdu/2*du1dx1];
        Jx2(:,:,i) = [ic' ix2' dfdx/2+dfdxdot/h+dfdu/2*du2dx2];
        Ju1(:,:,i) = [ic iu1 dfdu'/2];
        Ju2(:,:,i) = [ic iu2 dfdu'/2];
        JK1(:,:,i) = [ic' iKs1' dfdK/2+dfdu/2*du1dK1];
        JK2(:,:,i) = [ic' iKs2' dfdK/2+dfdu/2*du2dK2];
        
    end

    for i = 1:NperSU-1
        J(Jx1(:,1,i),Jx1(:,2,i)) = Jx1(:,3:4,i);
        J(Jx2(:,1,i),Jx2(:,2,i)) = Jx2(:,3:4,i);
        J(Ju1(:,1:2,i),Ju1(:,3,i)) = Ju1(:,4:5,i);
        J(Ju2(:,1:2,i),Ju2(:,3,i)) = Ju2(:,4:5,i);
        J(JK1(:,1,i),JK1(:,2,i)) = JK1(:,3:4,i);
        J(JK2(:,1,i),JK2(:,2,i)) = JK2(:,3:4,i);
    end
    
    if j == 1
        ix1 = (NperSU-1)*nvarpernode1+(1:nstates);
    else
        ix1 = nvarSU1+(j-2)*nvarSU+(NperSU-1)*nvarpernode+(1:nstates);
    end
    
    %Average of last nodes should equal pi/2, u should be zero
    J(end-1,ix1(1)) = 1/NSU;
    J(end,ix1(2)) = 1/NSU;
end
% 
% 
