function dobj = objgrad(X, params)

m = params.m;
l = params.l;
J = m*l^2;

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);

dobj = zeros(size(X));
for j = 1:params.NSU
    for i = 1:params.NperSU-1
        K = X(end-params.NperSU*2+i);
        Kd = X(end-params.NperSU+i);
        theta = X(ix(1));
        dtheta = X(ix(2));

        u = (X(iu)+K*theta+Kd*dtheta)/J;
        
        dobj(iu) = dobj(iu)+u/J;
        dobj(ix) = dobj(ix) + [u*K/J; u*Kd/J];
        dobj(end-params.NperSU*2+i) = dobj(end-params.NperSU*2+i)+u*theta/J;
        dobj(end-params.NperSU+i) = dobj(end-params.NperSU+i)+u*dtheta/J;
        if  theta > pi/2
            dobj(ix(1)) = dobj(ix(1))-0000*(pi/2-theta);
        end
        ix = ix+params.nvarpernode;
        iu = iu+params.nvarpernode;
    end

    %Final cost
    K = X(end-params.NperSU);
    Kd = X(end);
    theta = X(ix(1));
    dtheta = X(ix(2));
    u = (X(iu)+K*theta+Kd*dtheta)/J;
    
    dobj(ix) = dobj(ix) + [u*K/J; u*Kd/J];
    dobj(end-params.NperSU)= dobj(end-params.NperSU)+u*theta/J;
    dobj(end)= dobj(end)+u*dtheta/J;
    dobj(iu) = dobj(iu)+u/J;
    if  theta > pi/2
        dobj(ix(1)) = dobj(ix(1))-0000*(pi/2-theta);
    end
    
    ix = ix+params.nvarpernode;
    iu = iu+params.nvarpernode;
end

dobj = dobj/params.NSU;