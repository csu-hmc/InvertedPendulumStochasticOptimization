function [X, L, U] = getIniConBound(params, result)

% Bounds on x and u for each node
if nargin == 1
    Lpernode = [-2*pi;-1000;0;0;-1/2*params.m*params.g*params.l];
else
    Lpernode = [-2*pi;-1000;-inf;-inf;-1/2*params.m*params.g*params.l];
end
Upernode = [pi;1000;0;0;1/2*params.m*params.g*params.l];

LperSU = [];
UperSU = [];
for i = 1:params.NperSU
    LperSU = [LperSU;Lpernode];
    UperSU = [UperSU;Upernode];
end

L = [];
U = [];
for i = 1:params.NSU
    L = [L;LperSU];
    U = [U;UperSU];
end

if nargin == 1  
    %Mid initial guess
    X = 1/2*(L+U);
else
    X = [];
    for i = 1:params.NSU
        X = [X;result.X(1:result.params.nvarperSU)];
    end
end

    