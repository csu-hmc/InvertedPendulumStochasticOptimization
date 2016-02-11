function [X, L, U] = getIniConBound(params, result)

% Bounds on x and u for each node
Lpernode = [-3*pi;-1000;-100]; %-1/2*params.m*params.g*params.l
Upernode = [2*pi;1000;100]; %1/2*params.m*params.g*params.l

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
    L = [L;0;0]; %K and Kd
    U = [U;0;0];
else
    L = [L;0;0];%
    U = [U;-1000;-1000]; %-params.m*params.g*params.l;inf
end

if nargin == 1  
    %Mid initial guess
    X = 1/2*(L+U);
    X(end-1) = 0;% -params.m*params.g*params.l;
    X(end) = 0;%-10;
else
    X = [];
    for i = 1:params.NSU
        X = [X;result.X(1:result.params.nvarperSU)];
%         indices = 1+params.nvarpernode*(10:30);
%         X((i-1)*params.nvarperSU+indices) = -pi/2;
    end
    X = [X;result.X(end-1:end)];  
end

    