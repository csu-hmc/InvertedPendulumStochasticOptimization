function obj = objfun(X, params)

kT = params.ktheta;
kTd = params.kthetadot;
N = params.N;
m = params.m;
l = params.l;

J = m*l^2;

ix = 1:params.optstates;
iu = params.optstates+(1:params.ncontrols);

obj = 0;
for j = 1:params.NSU
    for i = 1:params.NperSU-1
        theta = X(ix(1));
        dtheta = X(ix(2));
        K = X(ix(3));
        Kd = X(ix(4));
        u = (X(iu)+K*theta+Kd*dtheta)/J;
        
        obj = obj+1/2*u^2;
        if  theta > (pi/2+0.01)
            obj = obj+1000*(pi/2+0.01-theta)^2;
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
    
    obj = obj + kT*(theta-pi/2)^2+kTd*dtheta^2+1/2*u^2;
    if theta > (pi/2+0.01)
        obj =  obj+1000*(pi/2+0.01-theta)^2;
    end
    ix = ix+params.nvarpernode;
    iu = iu+params.nvarpernode;
end

obj = obj/params.NSU;