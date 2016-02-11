function Ks = goodEigVals(params)

m = params.m;
l = params.l;
g = params.g;

J = m*l^2;

syms K Kd

lambda1 = (Kd/J+sqrt((Kd/J)^2+4*(m*g*l+K)/J))/2;
lambda2 = (Kd/J-sqrt((Kd/J)^2+4*(m*g*l+K)/J))/2;

eq1 = lambda1 <= 0;
eq2 = lambda2 <= 0;

solve([eq1, eq2], [K, Kd])