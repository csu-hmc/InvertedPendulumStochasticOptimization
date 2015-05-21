function [dobj1] = objgrad(X, params)

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

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);
N = params.N;

dobj1 = zeros(size(X));
for j = 1:params.NSU
    for i = 1:params.NperSU-1
        u = X(iu);
        
        dobj1(iu) = dobj1(iu)+u;
          if X(ix(1)) > (pi/2+0.1)
              dobj1(ix(1)) = dobj1(ix(1))-2000*(pi/2+0.1-X(ix(1)));
          end
        ix = ix+params.nvarpernode;
        iu = iu+params.nvarpernode;
    end
    %Final cost
    x = X(ix);
    u = X(iu);

    theta = x(1);
    dtheta = x(2);
    mhat = x(3);
    Sigma = x(4);
    
    dobj1(ix) = dobj1(ix) + [2*kT*mhat*(mhat*theta-pi/2)+2*theta*Sigma;
        2*kTd*dtheta;
        2*kT*theta*(mhat*theta-pi/2);
        theta^2];
    dobj1(iu) = dobj1(iu) + u;
end
