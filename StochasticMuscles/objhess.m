function H = objhess(X,sigma,lambda,params)

% keyboard;

itheta = params.ncontrols+params.ndof;
iact = params.ncontrols+params.ndof*2+(1:params.nmus);

H = spalloc(params.nvars,params.nvars,2*params.N*params.Nsamples);%length(X));
for j = 1:params.Nsamples
    for i = 1:params.N
%         theta = X(itheta);
%         u = X(iact); 
%         obj_eff = 1/2*sum(u.^2);
%         scalefun = 50*exp(u);
        H(iact,iact) = H(iact,iact)+eye(params.nmus);
        H(itheta,itheta) = H(itheta)+params.W*1;
        iact = iact+params.nvarpernode;
        itheta = itheta+params.nvarpernode;
    end
end
H = H;%/params.N/params.Nsamples*1e4;