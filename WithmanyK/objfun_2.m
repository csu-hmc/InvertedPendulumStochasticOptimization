function obj = objfun_2(X, params)

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);
iKs = params.nstates+params.ncontrols+(1:2);
nvarpernode1 = params.nvarpernode1;
nvarpernode = params.nvarpernode;

obj = 0;
for j = 1:params.NSU
    for i = 1:params.NperSU-1               
        x1 = X(ix);
		u1 = X(iu);
        K1 = X(iKs);
        u2 = X(iu+nvarpernode1);
        K2 = X(iKs+nvarpernode1);
        if j == 1
            x2 = X(ix+nvarpernode1);
        else
            x2 = X(ix+nvarpernode);
        end
%         
        x = (x1+x2)/2;
        u0 = (u1+u2)/2;
        Ks = (K1+K2)/2;
        u = findTorque(u0, Ks, x);

        obj = obj+1/2*u^2;

        if j == 1
            ix = ix+nvarpernode1;         
        else
            ix = ix+nvarpernode;
        end
        iu = iu+nvarpernode1;
        iKs = iKs+nvarpernode1;
    end
    
    % Restart iu and iKs
    iu = params.nstates+(1:params.ncontrols);
    iKs = params.nstates+params.ncontrols+(1:2);
end

obj = obj/params.NSU;