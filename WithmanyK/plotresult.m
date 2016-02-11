function [meanx, u0,K,Kd] = plotresult(result, plotyes)

if nargin == 1
    plotyes = 1;
end
params = result.params;
X = result.X;
X1 = X(1:params.nvarSU1);
x1 = reshape(X1,params.nvarpernode1,params.NperSU);
xperSU(:,:,1) = x1(1:2,:);
u0 = x1(3,:);
K = x1(4,:);
Kd = x1(5,:);
if params.NSU > 1
    Xelse = X(params.nvarSU1+1:end);
    xelse = reshape(Xelse,params.nvarpernode,params.NperSU,params.NSU-1);
    xperSU(:,:,2:params.NSU) = xelse;
end    

meanx = mean(xperSU,3);
stdx = std(xperSU,[],3);

uall = zeros(params.ncontrols,params.NperSU, params.NSU);
% Get complete input
for j = 1:params.NSU
    for i = 1:params.NperSU
        uall(:,i,j) = findTorque(u0(i), [K(i);Kd(i)],xperSU(:,i,j));
    end
end

uplot = zeros(params.NperSU,1);
for i = 1:params.NperSU
    uplot(i) = findTorque(u0(i), [K(i);Kd(i)],meanx(:,i));% u0(i)+[K(i) Kd(i)]*meanx(:,i);%(u0(i)+u0(i+1))/2;
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
    plot([0:h:T],uplot)
    hold on
    plot([0:h:T],uall(1,:,1),'r--')
    legend('Mean Torque', 'Swing Up 1')
    
    figure
    plot([0:h:T], K)
    hold on
    plot([0:h:T], Kd, 'r')
    legend('Position', 'Derivative')
end