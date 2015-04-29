function [dobj] = objgrad(X, params)

kT = params.ktheta;
kTd = params.kthetadot;

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);
N = params.N;

dobj = zeros(size(X));
for i = 1:N
    x = X(ix);
    u = X(iu);

    theta = x(1);
    dtheta = x(2);
    mhat = x(3);
    Sigma = x(4);

    dobj(ix) = dobj(ix)+ [2*kT*mhat*(mhat*theta-pi/2)+2*theta*Sigma;
        2*kTd*dtheta;
        2*kT*theta*(mhat*theta-pi/2);
        theta^2];
    dobj(iu) = dobj(iu)+u;
    
    ix = ix+params.nvarpernode;
    iu = iu+params.nvarpernode;
end

