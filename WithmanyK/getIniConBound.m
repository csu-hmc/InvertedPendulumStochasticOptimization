function [X, L, U] = getIniConBound(params, result)

% Bounds on x and u for each node
if params.NSU == 1 
    Lpernode1 = [-3*pi;-1000;-1000;0;0]; %theta thetadot u K Kd
    Upernode1 = [2*pi;1000;1000;0;0];
else
    Lpernode1 = [-3*pi;-1000;-1000;-1000;-1000]; %theta thetadot u K Kd
    Upernode1 = [2*pi;1000;1000;1000;1000];
end

Lpernode = [-3*pi;-1000];
Upernode = [2*pi;1000];

LperSU1 = [];
UperSU1 = [];
for i = 1:params.NperSU
    LperSU1 = [LperSU1;Lpernode1];
    UperSU1 = [UperSU1;Upernode1];
end

L = LperSU1;
U = UperSU1;
if params.NSU > 1
    LperSU = [];
    UperSU = [];
    for i = 1:params.NperSU
        LperSU = [LperSU;Lpernode];
        UperSU = [UperSU;Upernode];
    end
    for i = 1:params.NSU-1
        L = [L;LperSU];
        U = [U;UperSU];
    end
end

if nargin == 1  
    %Mid initial guess
    X = 1/2*(L+U);
else
    if result.params.nvars == params.nvars
        X = result.X;
    else
        oldX = reshape(result.X,params.nvarpernode1, params.NperSU);
        X = oldX(:);
        for i = 1:params.NSU-1
            addX = oldX(1:2,:); % Check if we use correct things
            X = [X;addX(:)];
        end
    end
end

% % Change this so it works with uu
% load('ResultInspect2.mat');
% params1 = result1.params;
% X = result1.X;
% X1 = X(1:params1.nvarSU1);
% x1 = reshape(X1,params1.nvarpernode1,params1.NperSU);
% xperSU(:,:,1) = x1(1:2,:);
% u0 = x1(3,:);
% K = x1(4,:);
% Kd = x1(5,:);
% if params1.NSU > 1
%     Xelse = X(params1.nvarSU1+1:end);
%     xelse = reshape(Xelse,params1.nvarpernode,params1.NperSU,params1.NSU-1);
%     xperSU(:,:,2:params1.NSU) = xelse;
% end    
% 
% 
% uall = zeros(params1.ncontrols,params1.NperSU, params1.NSU);
% % Get complete input
% for j = 1:params1.NSU
%     for i = 1:params1.NperSU-1
%         theta = (xperSU(1,i,j)+xperSU(1,i+1,j))/2;
%         dtheta = (xperSU(2,i,j)+xperSU(2,i+1,j))/2;
%         uall(:,i,j) = (u0(i)+u0(i+1))/2+(K(i)+K(i+1))/2*theta+(Kd(i)+Kd(i+1))/2*dtheta;
%     end
% end
% 
% X = [xperSU;uall;zeros(2,60)];
% X = X(:);