function x0 = findNeutralstate(params)

Lx = [0 ;-1000;0;0;0;0];
Ux = [pi; 1000;1;1;4;4];

x0 = fmincon(@obfun,zeros(params.nstates,1),[],[],[],[],Lx, Ux,@(x) neutfun(x,params));

end

function [ceq,dyns] = neutfun(x,params)
dyns = StocDyn(x, zeros(params.nstates,1), zeros(params.nmus,1), zeros(params.nmus,1), params);
ceq = [];

end

function obj = obfun(x)
obj = 1;
end