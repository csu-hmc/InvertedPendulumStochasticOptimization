function [meanx,u0, uall] = plotresult(result, plotyes)

if nargin == 1
    plotyes = 1;
end
params = result.params;
X = result.X;
X1 = X(1:params.nvarSU1);
x1 = reshape(X1,params.nvarpernode1,params.NperSU);
xperSU(:,:,1) = x1(1:params.nstates,:);
u0 = x1(params.nstates+(1:params.ncontrols),:);
if params.NSU > 1
    Xelse = X(params.nvarSU1+1:end-2);
    xelse = reshape(Xelse,params.nvarpernode,params.NperSU,params.NSU-1);
    xperSU(:,:,2:params.NSU) = xelse;
end

meanx = mean(xperSU,3);
stdx = std(xperSU,[],3);

ix = 1:params.nstates;
iu = params.nstates+(1:params.ncontrols);
uall = zeros(params.ncontrols,params.NperSU, params.NSU);
% Get complete input
K = X(end-1);
Kd = X(end);
for j = 1:params.NSU
    for i = 1:params.NperSU
        uall(:,i,j) = findTorque(u0(:,i),K,Kd,xperSU(1:params.ndof*2,i,j), params);
    end
end

for i = 1:60
    Fsee(:,i) = getMusDyns(meanx(:,i),zeros(6,1),uall(:,i),params);
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
    plot([0:h:T],u0(1,:))
    hold on
    plot([0:h:T],uall(1,:,1),'--')
    
    figure
    plot(Fsee')
    
    figure
    plot(meanx')
end