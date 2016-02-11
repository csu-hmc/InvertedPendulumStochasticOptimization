function [c] = ineqconfun(X, params,j)

NperSU = params.NperSU;
NSU = params.NSU;
nstates = params.nstates;
Xdes = params.Xdes;

x = reshape(X(1:end-2),params.nvarpernode, params.N);

xperSU = x(1:params.nstates,(j-1)*params.NperSU+(1:params.NperSU));

ix = 1:nstates/2;
ic = 1:nstates/2;

c = zeros(nstates/2*NperSU,1);
for i = 1:NperSU
    c(ic) = xperSU(1,i)-Xdes(i);
    ix = ix+params.nvarpernode;
    ic = ic+nstates/2;
end