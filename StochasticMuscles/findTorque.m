function [uout, dudx, dudK, dudKd] = findTorque(uin,Kval,Kdval,x, params)

K = sign(params.muscleparam.d')*Kval;
Kd= sign(params.muscleparam.d')*Kdval;

uout = uin+[K Kd]*(x-[pi/2;0]);

dudx = [K Kd];
dudK = (x(1)-pi/2)*sign(params.muscleparam.d');
dudKd= x(2)*sign(params.muscleparam.d');
