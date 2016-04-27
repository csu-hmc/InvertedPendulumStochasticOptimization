function obj = objfun(X, params)

% keyboard;

itheta = params.ncontrols+params.ndof;
iact = params.ncontrols+params.ndof*2+(1:params.nmus);

obj = 0;
for j = 1:params.Nsamples
    for i = 1:params.N
        theta = X(itheta);
        %Minimize muscle force
        u = X(iact); 
        obj_eff = 1/2*sum(u.^2);
%         scalefun = sum(50*exp(u));
        obj_track = 1/2*(theta-params.targetangle)^2;
        obj = obj+obj_eff+params.W*obj_track;
        iact = iact+params.nvarpernode;
        itheta = itheta+params.nvarpernode;
    end
end
%The mean activation
obj = obj/params.N/params.Nsamples;