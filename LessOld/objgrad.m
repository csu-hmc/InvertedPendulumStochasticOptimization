function dobj = objgrad(X, params)

kT = params.ktheta;
kTd = params.kthetadot;
m = params.m;
l = params.l;
J = m*l^2;

ix = 1:params.optstates;
iu = params.optstates+(1:params.ncontrols);

dobj = zeros(size(X));
for j = 1:params.NSU
    for i = 1:params.NperSU-1
        theta = X(ix(1));
        dtheta = X(ix(2));
        K = X(ix(3));
        Kd = X(ix(4));
        u = (X(iu)+K*theta+Kd*dtheta)/J;
        
        dobj(iu) = dobj(iu)+u/J;
        dobj(ix(1)) = dobj(ix(1))+u*K/J;
        dobj(ix(2)) = dobj(ix(2))+u*Kd/J;
        dobj(ix(3)) = dobj(ix(3))+u*theta/J;
        dobj(ix(4)) = dobj(ix(4))+u*dtheta/J;
        if  theta > (pi/2+0.01)
            dobj(ix(1)) = dobj(ix(1))-2000*(pi/2+0.01-theta);
        end
        ix = ix+params.nvarpernode;
        iu = iu+params.nvarpernode;
    end

    %Final cost
    theta = X(ix(1));
    dtheta = X(ix(2));
    K = X(ix(3));
    Kd = X(ix(4));
    u = (X(iu)+K*theta+Kd*dtheta)/J;
    
    dobj(ix) = dobj(ix) + [2*kT*(theta-pi/2)+u*K/J;
        2*kTd*dtheta+u*Kd/J;
        u*theta/J;
        u*dtheta/J];
        dobj(iu) = dobj(iu)+u/J;
    if  theta > (pi/2+0.01)
        dobj(ix(1)) = dobj(ix(1))-2000*(pi/2+0.01-theta);
    end
    
    ix = ix+params.nvarpernode;
    iu = iu+params.nvarpernode;
end

dobj = dobj/params.NSU;