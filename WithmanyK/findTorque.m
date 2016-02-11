function [uall, dudK, dudx] = findTorque(u0, K, x)

if nargin == 1
    uall = u0;
    dudK = [0 0];
    dudx = [0 0];
else
    uall = u0+K'*x;
    dudK = x';
    dudx = K';
end