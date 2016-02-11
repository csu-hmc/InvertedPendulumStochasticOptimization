function dobj = objgrad(X, params)

ncontrols = params.ncontrols;
nvarpernode1 = params.nvarpernode1;
nvarpernode = params.nvarpernode;
nstates = params.nstates;
nvarSU1 = params.nvarSU1;
nvarSU = params.nvarSU;

m = params.m;
l = params.l;
J = m*l^2;

dobj = zeros(size(X));
for j = 1:params.NSU
    for i = 1:params.NperSU
        if j == 1
            ix1 = (i-1)*nvarpernode1+(1:nstates);
        else
            ix1 = nvarSU1+(j-2)*nvarSU+(i-1)*nvarpernode+(1:nstates);
        end
        iu1 = (i-1)*nvarpernode1+nstates+(1:ncontrols);
        iKs1 = (i-1)*nvarpernode1+nstates+ncontrols+(1:2);
        
        x = X(ix1);
		u0 = X(iu1);
        Ks = X(iKs1);
        
        [u, dudK, dudx] = findTorque(u0, Ks, x);
        
        dobj(ix1) = dobj(ix1) + (u/J^2)*dudx';
        dobj(iu1) = dobj(iu1) + (u/J^2);
        dobj(iKs1)= dobj(iKs1)+ (u/J^2)*dudK';

    end

%     dobj(ix1) = dobj(ix1) + (u/J^2)*dudx';
%     dobj(iu1) = dobj(iu1) + (u/J^2);
%     dobj(iKs1)= dobj(iKs1)+ (u/J^2)*dudK';

end

dobj = dobj/params.NSU;