function [obj] = objfun(X, params)

ncontrols = params.ncontrols;
nvarpernode1 = params.nvarpernode1;
nvarpernode = params.nvarpernode;
nstates = params.nstates;
nvarSU1 = params.nvarSU1;
nvarSU = params.nvarSU;


m = params.m;
l = params.l;
J = m*l^2;

obj = 0;
for j = 1:params.NSU
    for i = 1:params.NperSU%-1  
        if j == 1
            ix = (i-1)*nvarpernode1+(1:nstates);
        else
            ix = nvarSU1+(j-2)*nvarSU+(i-1)*nvarpernode+(1:nstates);
        end
        iu = (i-1)*nvarpernode1+nstates+(1:ncontrols);
        iKs = (i-1)*nvarpernode1+nstates+ncontrols+(1:2);
        
        x = X(ix);
		u0 = X(iu);
        Ks = X(iKs);

        u = findTorque(u0, Ks, x)/J;

        obj = obj+1/2*u^2;
    end
end
obj = obj/params.NSU;