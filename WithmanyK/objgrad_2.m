function dobj = objgrad_2(X, params)

ix1 = 1:params.nstates;
iu1 = params.nstates+(1:params.ncontrols);
iKs1 = params.nstates+params.ncontrols+(1:2);
nvarpernode1 = params.nvarpernode1;
nvarpernode = params.nvarpernode;

dobj = zeros(size(X));
for j = 1:params.NSU
    for i = 1:params.NperSU-1               
        if j == 1
            ix2 = ix1 + nvarpernode1;         
        else
            ix2 = ix1 + nvarpernode;
        end
        iu2 = iu1 + nvarpernode1;
        iKs2= iKs1+ nvarpernode1;
        x1 = X(ix1);
		u1 = X(iu1);
		x2 = X(ix2);
		u2 = X(iu2);
        K1 = X(iKs1);
        K2 = X(iKs2);
        
        x = (x1+x2)/2;
        u0 = (u1+u2)/2;
        Ks = (K1+K2)/2;
        
        [u, dudK, dudx] = findTorque(u0, Ks, x);
        
        dobj(ix1) = dobj(ix1) + 1/2*u*dudx';%
        dobj(ix2) = dobj(ix2) + 1/2*u*dudx';
        dobj(iu1) = dobj(iu1) + 1/2*u;%
        dobj(iu2) = dobj(iu2) + 1/2*u;
        dobj(iKs1)= dobj(iKs1)+ 1/2*u*dudK';%
        dobj(iKs2)= dobj(iKs2)+ 1/2*u*dudK';

        if j == 1
            ix1 = ix1+nvarpernode1;            
        else
            ix1 = ix1+nvarpernode;
        end
        iu1 = iu1+nvarpernode1;
        iKs1 = iKs1+nvarpernode1;
    end
    
    % Restart iu and iKs
    iu1 = params.nstates+(1:params.ncontrols);
    iKs1 = params.nstates+params.ncontrols+(1:2);
end

dobj = dobj/params.NSU;