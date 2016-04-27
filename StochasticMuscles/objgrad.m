function dobj = objgrad(X, params)

% keyboard;

itheta = params.ncontrols+params.ndof;
iact = params.ncontrols+params.ndof*2+(1:params.nmus);

dobj = zeros(size(X));
for j = 1:params.Nsamples
    for i = 1:params.N
        theta = X(itheta);
        u = X(iact); 
        dobj(iact) = dobj(iact)+u;
        dobj(itheta) = dobj(itheta)+params.W*(theta-params.targetangle);
        iact = iact+params.nvarpernode;
        itheta = itheta+params.nvarpernode;
    end
end
dobj = dobj/params.N/params.Nsamples;