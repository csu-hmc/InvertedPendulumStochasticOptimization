function [J, params] = conjacstructure(L, U, params)
%---------------------------------------------------------------------------------
% determine sparsity structure of Jacobian

Jnnz = 1;						% temporary value needed for the first memory allocation
params.Jnnz = Jnnz;
nsame = 0;

if strcmp(params.optimizer, 'SNOPT')
    Jpattern = sparse(params.ncon+1,params.nvars);
    % do over dobj+jac
    % sparsity pattern needs to be the same 10 times in a row
    while (nsame < 10)
        X = L + (U-L).*randn(size(L));		% a random vector of unknowns
        dobj = objgrad(X, params);
        J = conjac(X, params);
        G = [dobj';J];
        newJpattern = double(G~=0) | Jpattern ;			% add any new nonzeros that were just found
        if (nnz(newJpattern - Jpattern) == 0)
            nsame = nsame + 1;
        else
            nsame = 0;
        end
        Jpattern = newJpattern ;
        Jnnz = nnz(Jpattern);
        params.Jnnz = Jnnz;
    end
else
    Jpattern = sparse(params.ncon,params.nvars);
    % sparsity pattern needs to be the same 10 times in a row
    while (nsame < 10)
        X = L + (U-L).*randn(size(L));		% a random vector of unknowns
        J = conjac(X, params);
        newJpattern = double(J~=0) | Jpattern ;			% add any new nonzeros that were just found
        if (nnz(newJpattern - Jpattern) == 0)
            nsame = nsame + 1;
        else
            nsame = 0;
        end
        Jpattern = newJpattern ;
        Jnnz = nnz(Jpattern);
        params.Jnnz = Jnnz;
    end
end

J = double(Jpattern);
params.Jnnz = Jnnz;
[params.DGrow,params.DGcol] = find(J ~= 0);