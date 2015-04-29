function [X, L, U] = getIniConBound(params)

% Bounds on x and u for each node
Lpernode = [-pi;-8;-2*pi;-100;-1/2*params.m*params.g*params.l];
Upernode = [pi;10;2*pi;100;1/2*params.m*params.g*params.l];

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

%Mid initial guess
X = 1/2*(L+U);

    