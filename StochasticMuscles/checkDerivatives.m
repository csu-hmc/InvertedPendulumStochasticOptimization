function [grad,grad_num,cjac, cjac_num] = checkDerivatives(X, params)

hh = 1e-6;

f = objfun(X, params);
grad = objgrad(X, params);
c = confun(X, params);
cjac = full(conjac(X, params));
cjac_num = zeros(params.ncon, params.nvars);
grad_num = zeros(params.nvars,1);
for i=1:params.nvars
    fprintf('checking objgrad and conjac for unknown %4d of %4d\n',i,params.nvars);
    Xisave = X(i);
    X(i) = X(i) + hh;
    cjac_num(:,i) = sparse(confun(X, params) - c)/hh;
    grad_num(i)   = (objfun(X, params) - f)/hh;
    hess_num(:,i) = sparse(objgrad(X, params) - grad)/hh;
    X(i) = Xisave;
end
		
% report maximal differences between analytical derivatives and numerical results
fprintf('Max. error in constraint jacobian: ');
matcompare(cjac, cjac_num);
fprintf('Max. error in objective gradient: ');
matcompare(grad, grad_num);

end

function matcompare(a,b);
	[maxerr,irow] = max(abs(a-b));
	[maxerr,icol] = max(maxerr);
	irow = irow(icol);
	fprintf('%9.6f at %d %d (%9.6f vs. %9.6f)\n', full(maxerr), irow, icol, full(a(irow,icol)), full(b(irow,icol)));
end