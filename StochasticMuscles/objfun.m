function obj = objfun(X, params)

itheta = params.ndof;
iact = params.ndof*2+(1:params.nmus);

obj = 0;
for j = 1:params.NSU
    for i = 1:params.NperSU
        theta = X(itheta);
        u = X(iact); 
        obj_eff = 1/2*sum(u.^2);
        obj_track = 1/2*(theta-pi/2)^2;
        obj = obj+obj_eff+params.W*obj_track;
        if j == 1
            iact = iact+params.nvarpernode1;
        else
            iact = iact+params.nvarpernode;
        end
    end
end

%The mean activation
% obj = obj/params.NSU/params.NperSU;