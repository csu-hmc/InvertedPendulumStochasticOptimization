function [obj,obj_eff_tot,obj_track_tot] = objfun(X, params)

% keyboard;

itheta = params.ncontrols+params.ndof;
iact = params.ncontrols+params.ndof*2+(1:params.nmus);

obj = 0;
obj_track_tot = 0;
obj_eff_tot=0;
for j = 1:params.Nsamples
    for i = 1:params.N
        theta = X(itheta);
        
        u = findTorque(X(1:params.ncontrols),X(end-params.Ks+1:end),X(itheta:itheta+1),params);
        a = X(iact);
        
        obj_eff = 1/2*sum(u.^2);
%         obj_eff = 1/2*sum(a.^2);
        obj_track = 1/2*(theta-params.targetangle)^2*1e3;
        obj = obj+obj_eff+params.W*obj_track;%obj+obj_track;%

        obj_eff_tot = obj_eff_tot+obj_eff;
        obj_track_tot = obj_track_tot+obj_track;
        itheta = itheta+params.nvarpernode;
    end
end
%The mean activation
obj = obj/params.N/params.Nsamples;