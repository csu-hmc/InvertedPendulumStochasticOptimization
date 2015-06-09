function obj = objfun(X, params)

m = params.m;
l = params.l;
J = m*l^2;

K = X(end-1);
Kd = X(end);

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);

obj = 0;
for j = 1:params.NSU
    for i = 1:params.NperSU-1
        theta = X(ix(1));
        dtheta = X(ix(2));
        u = (X(iu)+K*theta+Kd*dtheta)/J;
        
        obj = obj+1/2*u^2;
        if  theta > pi/2
            obj = obj+0000*(pi/2-theta)^2;
        end
        ix = ix+params.nvarpernode;
        iu = iu+params.nvarpernode;
    end
    %Final cost
    theta = X(ix(1));
    dtheta = X(ix(2));
    u = (X(iu)+K*theta+Kd*dtheta)/J;
    
    obj = obj + 1/2*u^2;
    if theta > pi/2
        obj =  obj+0000*(pi/2-theta)^2;
    end
    ix = ix+params.nvarpernode;
    iu = iu+params.nvarpernode;
end

obj = obj/params.NSU;