clear all
% close all
clc

xx = [1 2 3 4 5 6 7 8 9 10];% 20 30 40 50 60 70 80 90];% 50 60 70 80 90 100 110 130 150 170 190 200 300 400 500 600 700 800];

RMS = zeros(length(xx),10);
obj = zeros(length(xx),10);
K = zeros(2,length(xx),10);
meanx = zeros(length(xx),60,10);
meanxdot = zeros(length(xx),60,10);
for i = 1:length(xx)
    for j = 1:10
        loadname = ['Solutions/WithConstraint/Solution_' num2str(xx(i)) 'SUs' num2str(j) '.mat'];
        load(loadname)
        obj(i,j) = result.obj;
        K(:,i,j) = result.X(end-1:end);
        [uperSU, uall,meanxx] = plotresult(result, 0);
        meanx(i,:,j) = meanxx(1,:);
        meanxdot(i,:,j) = meanxx(2,:);
        rmsu = rms(uall);
        RMS(i,j) = mean(rmsu(:));
    end    
end

K1 = mean(K,3);
K1std = std(K,[],3);
obj1 = mean(obj,2);
RMS1 = mean(RMS,2);
figure
errorbar(xx,K1(1,:), K1std(1,:))
hold on
errorbar(xx,K1(2,:), K1std(2,:), 'r')
xlabel('Number of Swing-Ups')
ylabel('Feedback Gain')
legend('Kp', 'Kd')
figure
errorbar(xx,obj1, std(obj,[],2))
xlabel('Number of Swing-Ups')
ylabel('Objective of Optimal Solution')
figure
errorbar(xx,RMS1, std(RMS,[],2))
T = result.params.T;
h = T/(result.params.NperSU-1);
figure
xdots = mean(meanxdot,3);
plot(xx,xdots(:,end))
xlabel('Number of Swing-Ups')
ylabel('Average Final Velocity')
