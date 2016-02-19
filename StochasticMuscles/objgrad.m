function dobj = objgrad(X, params)

itheta = params.ndof;
iact = params.ndof*2+(1:params.nmus);

dobj = zeros(size(X));
for j = 1:params.NSU
    for i = 1:params.NperSU
        theta = X(itheta);
        u = X(iact); 
        dobj(iact) = dobj(iact)+u;
        dobj(itheta) = dobj(itheta)+params.W*(theta-pi/2);
        if j == 1
            iact = iact+params.nvarpernode1;
        else
            iact = iact+params.nvarpernode;
        end
    end
end

% dobj = dobj/params.NSU/params.NperSU;