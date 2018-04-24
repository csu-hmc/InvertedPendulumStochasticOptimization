function [X, L, U] = getIniConBound(params, nou0, result)

if nargin == 1
    nou0 = 0;
end

% Bounds on state

% if params.asat
%     Lx = [-pi;-1000;-10;-10;0;0];
%     Ux = [2*pi; 1000;100;100;4;4];
% else
    Ux = [6*pi; 1000;2;2;1.5;1.5];
    Lx = [-5*pi;-1000;0;0;0.5;0.5];
% end
X0x = params.xneutral;

% Bounds on controls, with or without co-contraction 
if nou0 == 1
    Lu = 1e-4+zeros(params.ncontrols,1);
    Uu = repmat(1e-2,params.ncontrols,1);
else
    Lu = repmat(0,params.ncontrols,1);
    Uu = repmat(0.4,params.ncontrols,1);
end
X0u = 1e-4+zeros(params.ncontrols,1);

L = [Lu; repmat(Lx,params.N*params.Nsamples,1)];
U = [Uu; repmat(Ux,params.N*params.Nsamples,1)];
X0 = [X0u;repmat(X0x,params.N*params.Nsamples,1)];

if params.omega == 0
    L = [L;zeros(params.Ks,1)]; %K and Kd
    U = [U;zeros(params.Ks,1)];
    %Mid initial guess
    X = [X0;zeros(params.Ks,1)];
%     U(1:params.ncontrols)  = ones(params.ncontrols,1);
else
    X = [X0;zeros(params.Ks,1)];
    L = [L;repmat([-10;0],params.Ks/2,1)];
    U = [U;repmat([0;10],params.Ks/2,1)];
end

if nargin == 3
    if result.params.N == params.N
        if params.Nsamples ~= result.params.Nsamples
            X = [result.X(1:params.ncontrols);repmat(result.X(params.ncontrols+(1:result.params.nvarSU1)),params.Nsamples,1);result.X(end-params.Ks+1:end)];
        else
            X = result.X;
        end
    else
        X(1:params.ncontrols) = result.X(1:params.ncontrols);
        X(end-params.Ks+1:end) = result.X(end-params.Ks+1:end);
    end
    %hacks
    X(1:2) = [0.01; 0.01];%rand(2,1);%ones(2,1);% [0.5;0.5];%
    X(end-params.Ks+1:end) = [-0.01; .01;-0.01;0.01];%.*rand(4,1);
end
