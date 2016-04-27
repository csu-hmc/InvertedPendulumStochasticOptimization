function [f,G] = objconfunp(X)

global params

obj = objfun(X, params);
dobj = objgrad(X, params);

c = confun(X, params);
J = conjac(X, params);

f = [obj;c];
G = [dobj';J];
% norm(c)
% norm(f)