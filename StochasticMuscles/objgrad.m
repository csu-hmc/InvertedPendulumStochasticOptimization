function dobj = objgrad(X, params)

% keyboard;
ix = params.ncontrols+(1:params.ndof*2);
iact = params.ncontrols+params.ndof*2+(1:params.nmus);

dobj = zeros(size(X));
for j = 1:params.Nsamples
    for i = 1:params.N
        theta = X(ix(1));
        
        [u, dudx, dudK, dudKd] = findTorque(X(1:params.ncontrols),X(end-params.Ks+1:end),X(ix),params);
%         a = X(iact);
         
        % Deviation
        dobj(ix(1)) = dobj(ix(1))+params.W*(theta-params.targetangle)*1e3;
        
        % Muscles
%         dobj(iact) = dobj(iact)+a;
        dobj(ix) = dobj(ix)+dudx*u;
        dobj(end-params.Ks+1:end) = dobj(end-params.Ks+1:end)+[dudK*u;dudKd*u];
        dobj(1:params.ncontrols) = dobj(1:params.ncontrols)+u;
        ix = ix+params.nvarpernode;
    end
end
dobj = dobj/params.N/params.Nsamples;