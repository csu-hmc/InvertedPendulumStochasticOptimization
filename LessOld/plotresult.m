function plotresult(result, params)

X = result.X;
x = reshape(X(1:end),params.nvarpernode, params.N);

xperSU = zeros(params.nstates, params.NperSU, params.NSU);
uperSU = zeros(params.ncontrols,params.NperSU, params.NSU);
for i = 1:params.NSU
    xperSU(:,:,i) = x(1:params.nstates,(i-1)*params.NperSU+(1:params.NperSU));
    uperSU(:,:,i) = x(params.optstates+1:end,(i-1)*params.NperSU+(1:params.NperSU));
end

meanx = mean(xperSU,3);
meanu = mean(uperSU,3);

figure
plot(meanx(1,:));
hold on
for i = 1:size(xperSU,3)
    plot(xperSU(1,:,i),'color',0.5*[1 1 1])
end
plot([0 params.NperSU], [pi/2+0.01 pi/2+0.01], 'r')
figure
plot(meanu(1,:))