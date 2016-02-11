KK = 10*rand(2,60);

params = result1.params;
X = result1.X;
X1 = X(1:params.nvarSU1);
x1 = reshape(X1,params.nvarpernode1,params.NperSU);
xperSU(:,:,1) = x1(1:2,:);
if params.NSU > 1
    Xelse = X(params.nvarSU1+1:end);
    xelse = reshape(Xelse,params.nvarpernode,params.NperSU,params.NSU-1);
    xperSU(:,:,2:params.NSU) = xelse;
end    

uall = x1(3,:);
K = KK(1,:);
Kd = KK(2,:);
u0 = zeros(params.ncontrols,params.NperSU);
fu = zeros(params.ncontrols,params.NperSU);
% Get complete input
for i = 1:params.NperSU
    theta = xperSU(1,i,1);
    dtheta = xperSU(2,i,1);
    fu(:,i) = K(i)*theta+Kd(i)*dtheta;
    u0(:,i) = uall(i)-fu(:,i);
end

xwithK = [x1(1:2,:);u0;K;Kd];
XwithK = xwithK(:);
% fu1 = fu/params.m/params.l^2;
[cons2,u2] = confun(result1.X,result1.params);
[cons1,u1] = confun(XwithK,result1.params);
% u01 = u0/params.m/params.l^2;
% uall1 = uall/params.m/params.l^2;