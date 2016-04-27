function [uout, duoutdx, duoutdu, duoutdK, duoutdKd] = findTorque(uin,Ks,x, params)

if params.ncontrols == 1
    K = sign(params.muscleparam.d')*Ks(1);
    Kd= sign(params.muscleparam.d')*Ks(2);
else
    K = Ks(1:2);
    Kd= Ks(3:4);
end

u = uin+[K Kd]*10*(x-[params.targetangle;0]);
if params.asat
    upos = 1/2*(u+sqrt(u.^2+params.epsilon^2));
    uout = 1-1/2*((1-upos)+sqrt((1-upos).^2+params.epsilon^2));
    duposdu = diag(1/2+1/2*(u.^2+params.epsilon^2).^(-1/2).*u);
    duoutdupos = diag(1/2+1/2*((1-upos).^2+params.epsilon^2).^(-1/2).*(1-upos));
    duoutdu = duoutdupos*duposdu;
else
    uout = u;
    duoutdu = eye(2);
end
duoutdx = duoutdu*10*[K Kd];
if params.ncontrols == 1
    duoutdK = duoutdu*10*((x(1)-params.targetangle)*sign(params.muscleparam.d'));
    duoutdKd= duoutdu*10*(x(2)*sign(params.muscleparam.d'));
else
    duoutdK = duoutdu*10*(x(1)-params.targetangle);
    duoutdKd= duoutdu*10*x(2);
end