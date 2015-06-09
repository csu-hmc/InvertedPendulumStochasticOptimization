function [uperSU, uall,meanx] = plotresult(result, plotyes)

if nargin == 1
    plotyes = 1;
end
params = result.params;
X = result.X;
x = reshape(X(1:end-2),params.nvarpernode, params.N);

xperSU = zeros(params.nstates, params.NperSU, params.NSU);
uperSU = zeros(params.ncontrols,params.NperSU, params.NSU);
for i = 1:params.NSU
    xperSU(:,:,i) = x(1:params.nstates,(i-1)*params.NperSU+(1:params.NperSU));
    uperSU(:,:,i) = x(params.nstates+1:end,(i-1)*params.NperSU+(1:params.NperSU));
end

meanx = mean(xperSU,3);
meanu = mean(uperSU,3);

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);
uall = zeros(params.ncontrols,params.NperSU, params.NSU);
% Get complete input
K = X(end-1);
Kd = X(end);
for j = 1:params.NSU
    for i = 1:params.NperSU
        theta = X(ix(1));
        dtheta = X(ix(2));
        uall(:,i,j) = X(iu)+K*theta+Kd*dtheta;

        ix = ix+params.nvarpernode;
        iu = iu+params.nvarpernode;
    end
end

if plotyes == 1   
    figure
    hold on
    for i = 1:size(xperSU,3)
        plot(xperSU(1,:,i),'color',0.5*[1 1 1])
    end
    plot(meanx(1,:), 'b', 'LineWidth', 1.5);
    plot([0 params.NperSU], [pi/2 pi/2], 'r')
    figure
    plot(meanu(1,:))

    hold on
    plot(uall(1,:,1),'--')
end