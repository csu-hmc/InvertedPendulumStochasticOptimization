function [X, L, U] = getIniConBound(params, result)

% Bounds on state
Lx = [1/9*pi;-1000;0;0;0.5;0.5];
Ux = [8/9*pi; 1000;2;2;1.5;1.5];

% Bounds on controls
Lu = [-100; -100];
Uu = [ 101;  101];

% Bounds on x and u for each node
Lpernode = [Lx;Lu];
Upernode = [Ux;Uu];

LperSU = repmat(Lpernode,params.NperSU,1);
UperSU = repmat(Upernode,params.NperSU,1);

L = LperSU;
U = UperSU;
if params.NSU > 1
    LperSU = [];
    UperSU = [];
    for i = 1:params.NperSU
        LperSU = [LperSU;Lx];
        UperSU = [UperSU;Ux];
    end
    for j = 2:params.NSU
        L = [L;LperSU];
        U = [U;UperSU];
    end
end

if nargin == 1  
    L = [L;0;0]; %K and Kd
    U = [U;0;0];
    %Mid initial guess
    X = 1/2*(L+U);
else
    if result.params.nvars == params.nvars
        X = result.X;
    else
        oldX = reshape(result.X(1:end-2),params.nvarpernode1, params.NperSU);
        X = oldX(:);
        for i = 1:params.NSU-1
            addX = oldX(1:params.nstates,:); % Check if we use correct things
            X = [X;addX(:)];
        end
        X = [X;result.X(end-1:end)];
    end
    L = [L;-1000;-1000];%
    U = [U;0;0]; %-params.m*params.g*params.l;inf
end

    