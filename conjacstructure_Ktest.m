function [J, params] = conjacstructure_Ktest(L, U, params)
%---------------------------------------------------------------------------------
% determine sparsity structure of Jacobian

Jnnz = 1;						% temporary value needed for the first memory allocation
params.Jnnz = Jnnz;
Jpattern = sparse(params.ncon,params.nvars);
nsame = 0;
% sparsity pattern needs to be the same 10 times in a row
while (nsame < 10)
    X = L + (U-L).*randn(size(L));		% a random vector of unknowns
    J = conjac_Ktest(X, params);
    newJpattern = double(J~=0) | Jpattern ;			% add any new nonzeros that were just found
    if (nnz(newJpattern - Jpattern) == 0)
        nsame = nsame + 1;
    else
        nsame = 0;
    end
    Jpattern = newJpattern ;
    Jnnz = nnz(Jpattern);
    params.Jnnz = Jnnz;
    params.Jpattern = Jpattern;
end

J = double(Jpattern);
params.Jnnz = Jnnz;