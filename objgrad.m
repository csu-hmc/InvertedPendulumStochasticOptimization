function dobj = objgrad(X, params)

m = params.m;
l = params.l;
J = m*l^2;

K = X(end-1);
Kd = X(end);

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);

dobj = zeros(size(X));
for j = 1:params.NSU
    for i = 1:params.NperSU-1
        theta = X(ix(1));
        dtheta = X(ix(2));

        u = (X(iu)+K*theta+Kd*dtheta)/J;
        
        dobj(iu) = dobj(iu)+u/J;
        dobj(ix) = dobj(ix) + [u*K/J; u*Kd/J];
        dobj(end-1) = dobj(end-1)+u*theta/J;
        dobj(end) = dobj(end)+u*dtheta/J;
        if  theta > pi/2
            dobj(ix(1)) = dobj(ix(1))-0000*(pi/2-theta);
        end
        ix = ix+params.nvarpernode;
        iu = iu+params.nvarpernode;
    end

    %Final cost
    theta = X(ix(1));
    dtheta = X(ix(2));
    u = (X(iu)+K*theta+Kd*dtheta)/J;
    
    dobj(ix) = dobj(ix) + [u*K/J; u*Kd/J];
    dobj(end-1)= dobj(end-1)+u*theta/J;
    dobj(end)= dobj(end)+u*dtheta/J;
    dobj(iu) = dobj(iu)+u/J;
    if  theta > pi/2
        dobj(ix(1)) = dobj(ix(1))-0000*(pi/2-theta);
    end
    
    ix = ix+params.nvarpernode;
    iu = iu+params.nvarpernode;
end

dobj = dobj/params.NSU;