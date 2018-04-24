function [u, dudx, dudK, dudKd] = findTorque(uin,Ks,x, params)

if params.ncontrols == 1
    K = sign(params.muscleparam.d')*Ks(1);
    Kd= sign(params.muscleparam.d')*Ks(2);
else
    K = Ks(1:2)*100;
    Kd= Ks(3:4)*10;
end

u = uin+[K Kd]*(x-[params.targetangle;0]);

dudx = [K'; Kd'];
if params.ncontrols == 1
    dudK = diag(((x(1)-params.targetangle)*sign(params.muscleparam.d')));
    dudKd= diag((x(2)*sign(params.muscleparam.d')));
else
    dudK = diag((x(1)-params.targetangle)*[1;1]*100);
    dudKd= diag(x(2)*[1;1]*10);
end