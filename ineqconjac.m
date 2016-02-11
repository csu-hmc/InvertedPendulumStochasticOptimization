function J = ineqconjac(X, params,j)

NperSU = params.NperSU;
NSU = params.NSU;
nstates = params.nstates;
nvars = params.nvars;

ix = params.nconperSU*(j-1)+(1:nstates/2);
ic = 1:nstates/2;

J = spalloc(nstates/2*NperSU,nvars,params.N);

ix1 = 1:nstates/2;
for i = 1:NperSU
    J(ic,ix) = 1;
    ix2 = ix1;
    for j = 1:params.NSU
        J(ic,ix2) = J(ic,ix2)-1/params.NSU;
        ix2 = ix2+params.nvarperSU;
    end
    ix1 = ix1+params.nvarpernode;
    ix = ix+params.nvarpernode;
    ic = ic+nstates/2;
end