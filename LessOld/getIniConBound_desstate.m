function [X, L, U] = getIniConBound(params, result)

% Bounds on x and u for each node
Lpernode = [-pi;-10;-6;-100;-1/2*params.m*params.g*params.l;-1/2*params.m*params.g*params.l];
Upernode = [pi;10;6;100;1/2*params.m*params.g*params.l;1/2*params.m*params.g*params.l];

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

% Add feedback gain to state
L = [L;-100;-10];
U = [U;0;0];

if nargin == 1  
    %Mid initial guess
    X = 1/2*(L+U);
else
    X = [];
    for i = 1:params.NSU
        X = [X;result.X(1:result.params.nvarperSU)];
    end
    X = [X;-5;-1];
end

    