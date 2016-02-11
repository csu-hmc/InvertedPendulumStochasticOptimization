clear all
close all
clc

xx = [1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100];

RMSu = zeros(length(xx),10);
obj = zeros(length(xx),10);
KK = zeros(120,length(xx),10);
meanx = zeros(length(xx),60,10);
meanxdot = zeros(length(xx),60,10);
for i = 1:length(xx)
    for j = 1:10
        loadname = ['Solutions/Solution_' num2str(xx(i)) 'SUs' num2str(j) '.mat'];
        load(loadname)
        obj(i,j) = result(end).obj;
%         K(:,i,j) = result(5).X(end-1:end);
        [xnow, u0,K,Kd] = plotresult(result(end), 0);
        for k = 1:length(u0)
            uall(k) = findTorque(u0(k), [K(k);Kd(k)], xnow(:,k));
        end
        KK(:,i,j) = [K Kd];
        meanx(i,:,j) = xnow(1,:);
        meanxdot(i,:,j) = xnow(2,:);
        rmsu = rms(uall);
        uRMS(i,j) = mean(rmsu(:));
    end    
end

% K1 = mean(K,3);
% K1std = std(K,[],3);
obj1 = mean(obj,2);
RMS1 = mean(uRMS,2);
% figure
% errorbar(xx,K1(1,:), K1std(1,:))
% hold on
% errorbar(xx,K1(2,:), K1std(2,:), 'r')
% xlabel('Number of Swing-Ups')
% ylabel('Feedback Gain')
% legend('Kp', 'Kd')
figure
errorbar(xx,obj1, std(obj,[],2))
xlabel('Number of Swing-Ups')
ylabel('Objective of Optimal Solution')
figure
errorbar(xx,RMS1, std(uRMS,[],2))
T = result(end).params.T;
h = T/(result(end).params.NperSU-1);
% figure
% xdots = mean(meanxdot,3);
% plot(xx,xdots(:,end))
% xlabel('Number of Swing-Ups')
% ylabel('Average Final Velocity')

%% Check increasing noise

RMSu = zeros(length(result),10);
obj = zeros(length(result),10);
K = zeros(2,length(result),10);
meanx = zeros(length(result),60,10);

for j = 1:10
    loadname = ['Solutions/Solution_100SUs' num2str(j) '.mat'];
    load(loadname)
    for i = 1:length(result)
        obj(i,j) = result(i).obj;
        [xnow, u0,K,Kd] = plotresult(result(i), 0);
        for k = 1:length(u0)
            uall(k) = findTorque(u0(k), [K(k);Kd(k)], xnow(:,k));
        end
        meanx(i,:,j) = xnow(1,:);
        rmsu = rms(uall);
        RMSu(i,j) = mean(rmsu(:));
    end
end
figure
meantrajs = mean(meanx,3);
stdtrajs = std(meanx,[],3);
plotstds = meantrajs-stdtrajs;
plotstds(:,end+1:2*end) = meantrajs(:,end:-1:1)+stdtrajs(:,end:-1:1);
xstds = [0:h:T T:-h:0];
plot([0:h:T], meantrajs', 'linewidth', 1.5)
legend('\sigma = 0', '\sigma = 1e-2', '\sigma = 1e-1', '\sigma = 1', 'Location', 'NorthWest', 'FontSize', 12)
xlabel('Time [s]')
ylabel('Pendulum Angle [rad]')
figure
fill(xstds,plotstds(2,:),[1 0.9 0.9])
hold on
fill(xstds,plotstds(3,:),[0.9 0.9 1])
fill(xstds,plotstds(4,:),[0.9 1 0.9])
legend('\sigma = 1e-3', '\sigma = 1e-1', '\sigma = 1', 'Location', 'NorthWest')
xlabel('Time [s]')
ylabel('Pendulum Angle [\theta]')
