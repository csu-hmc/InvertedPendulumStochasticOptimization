function [obj1] = objfun(X, params)

kT = params.ktheta;
kTd = params.kthetadot;
N = params.N;

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);

obj = 0;
for i = 1:N
    x = X(ix);
    u = X(iu);

    theta = x(1);
    dtheta = x(2);
    mhat = x(3);
    Sigma = x(4);

    obj = obj + kT*((mhat*theta-pi/2)^2+theta^2*Sigma)+kTd*dtheta^2+1/2*u^2;
    
    ix = ix+params.nvarpernode;
    iu = iu+params.nvarpernode;
end

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);
% 
obj1 = 0;
for j = 1:params.NSU
    for i = 1:params.NperSU-1
        u = X(iu);
        
        obj1 = obj1+1/2*u^2;
          if X(ix(1)) > (pi/2+0.1)
              obj1 = obj1+1000*(pi/2+0.1-X(ix(1)))^2;
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
    
    obj1 = obj1 + kT*((mhat*theta-pi/2)^2+theta^2*Sigma)+kTd*dtheta^2+1/2*u^2;
end