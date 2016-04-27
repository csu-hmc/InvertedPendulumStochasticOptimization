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
icineq = params.nconeq+(1:params.nmus);
        
nvarpernode = params.nvarpernode;
ncon = params.ncon;
ndof = params.ndof;

c = zeros(ncon,1);

for j = 1:params.Nsamples
    % % First node should match initial conditions
    x1 = X(ix);
    c(ic) = x1-params.xneutral;% [x1(1)-pi/2; x1(2)];%;  x1(5:6)-1];%x1(3:4)-u;
    ic = ic+nstates;%+ndof*2; 

    % Constraints on dynamics
    for i = 1:N-1
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
        %c(end)+x1(1)/N/params.Nsamples;
%         c(icineq) = u1;

        ix = ix+nvarpernode;
        ic = ic+nstates;
%         icineq = icineq+params.nmus;
    end

    x1 = X(ix);
    u1 = findTorque(u,X(end-params.Ks+1:end),x1(1:ndof*2),params);

%     c(icineq) = u1;
%     c(end) = c(end)+x1(1)/params.Nsamples/N; %
%     icineq = icineq+params.nmus;
    ix = ix + nvarpernode;
end
c(end) = (u(1)-u(2))/u(1);% c(end) = c(end)-params.targetangle;