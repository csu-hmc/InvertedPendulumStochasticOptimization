function [c] = confun(X, params)
% keyboard

% define constraints

nstates = params.nstates;
ncontrols = params.ncontrols;
h = params.h;
N = params.N;

ix = ncontrols+(1:nstates);
u = X(1:ncontrols); %same open-loop signal
ic = 1:nstates;
        
nvarpernode = params.nvarpernode;
ncon = params.ncon;
ndof = params.ndof;

c = zeros(ncon,1);

for j = 1:params.Nsamples
    % Constraints on dynamics
    for i = 1:N
        x1 = X(ix);
        if i == N
            x2 = X(ncontrols+(1:nstates));
        else
            x2 = X(ix+nvarpernode);
        end

        u1 = findTorque(u,X(end-params.Ks+1:end),x1(1:ndof*2),params);
        u2 = findTorque(u,X(end-params.Ks+1:end),x2(1:ndof*2),params);

        omega = params.omega(:,(j-1)*params.N+i);
        if strcmp(params.method, 'ME')
            dyns = StocDyn((x1+x2)/2,(x2-x1)/h,(u1+u2)/2, omega, params);
        else
            dyns = StocDyn(x2,(x2-x1)/h,u2, omega, params);
        end
        c(ic) = dyns;

        ix = ix+nvarpernode;
        ic = ic+nstates;
    end
end
% c(end) = (u(1)-u(2))/u(1);%