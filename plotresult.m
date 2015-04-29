function plotresult(result, params)

X = result.X;
x = reshape(X,params.nvarpernode, params.N);

xperSU = zeros(params.nstates, params.NperSU, params.NSU);
uperSU = zeros(1,params.NperSU, params.NSU);
for i = 1:params.NSU
    xperSU(:,:,i) = x(1:params.nstates,(i-1)*params.NperSU+(1:params.NperSU));
    uperSU(:,:,i) = x(end,(i-1)*params.NperSU+(1:params.NperSU));
end

meanx = mean(xperSU,3);
meanu = mean(uperSU,3);

plot(meanx(1,:));
figure
plot(meanu)
figure
plot(meanx(4,:));
figure
plot(meanx(3,:));