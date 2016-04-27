function [H, params] = hessstructure(L, U, params)
%---------------------------------------------------------------------------------
% determine sparsity structure of Jacobian

Hnnz = 1;						% temporary value needed for the first memory allocation
params.Hnnz = Hnnz;
Hpattern = sparse(params.nvars,params.nvars);
nsame = 0;
% sparsity pattern needs to be the same 10 times in a row
while (nsame < 3)
    X = L + (U-L).*randn(size(L));		% a random vector of unknowns
    H = objhess(params);
    newHpattern = double(H~=0) | Hpattern ;			% add any new nonzeros that were just found
    if (nnz(newHpattern - Hpattern) == 0)
        nsame = nsame + 1;
    else
        nsame = 0;
    end
    Hpattern = newHpattern ;
    Hnnz = nnz(Hpattern);
    params.Hnnz = Hnnz;
    params.Hpattern = Hpattern;
end

H = double(Hpattern);
params.Hnnz = Hnnz;