function [meanu, uall,meanx] = plotresult(result, plotyes)

if nargin == 1
    plotyes = 1;
end
params = result.params;
X = result.X;
x = reshape(X(1:end-params.Ks),params.nvarpernode, params.N);

xperSU = zeros(params.nstates, params.NperSU, params.NSU);
uperSU = zeros(params.ncontrols,params.NperSU, params.NSU);
for i = 1:params.NSU
    xperSU(:,:,i) = x(1:params.nstates,(i-1)*params.NperSU+(1:params.NperSU));
    uperSU(:,:,i) = x(params.nstates+1:end,(i-1)*params.NperSU+(1:params.NperSU));
end

meanx = mean(xperSU,3);
stdx = std(xperSU,[],3);
meanu = mean(uperSU,3);

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);
uall = zeros(params.ncontrols,params.NperSU, params.NSU);
% Get complete input
if isfield(params, 'Ks')
    K = X(end-119:end-60);
    Kd = X(end-59:end);
else
    K = X(end-1)+zeros(60,1);
    Kd = X(end)+zeros(60,1);
end
for j = 1:params.NSU
    for i = 1:params.NperSU
        theta = X(ix(1));
        dtheta = X(ix(2));
        uall(:,i,j) = X(iu)+K(i)*theta+Kd(i)*dtheta;

        ix = ix+params.nvarpernode;
        iu = iu+params.nvarpernode;
    end
end

if plotyes == 1
    T = result.params.T;
    h = T/(result.params.NperSU-1);
    xstds = [0:h:T T:-h:0];
    plotstds = meanx-stdx;
    plotstds(:,end+1:2*end) = meanx(:,end:-1:1)+stdx(:,end:-1:1);
    figure
%     fill(xstds,plotstds(1,:),[1 0.9 0.9])
    hold on
    for i = 1:size(xperSU,3)
        plot([0:h:T],xperSU(1,:,i),'color',0.2*[1 1 1])
    end
    plot([0:h:T],meanx(1,:), 'b', 'LineWidth', 1.5);
    plot([0 T], [pi/2 pi/2], 'r')
    
    figure
    plot([0:h:T],meanu(1,:))
    hold on
    plot([0:h:T],uall(1,:,1),'--')
end