clear all
close all
clc

% Weights = [0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100];
% Times = [1 3 5 8 10 25 50 75 100 250 500 750];% 1000];
Ns = [100 200 300 400 500 600 700 800 900 1000];%To use for general parameters.
% Ts = [1 5 10 20 30 40 50];% 100 200];

u0 = zeros(2,length(Ns),5);
Ks = zeros(2,length(Ns),5);
obj = zeros(length(Ns),5);
for i = 1:length(Ns)
    for j = 1:5
        
        filename = ['Solutions/Results_', num2str(Ns(i)), 'hcheck10_' num2str(j) '.mat'];
%         filename = ['Solutions/Results_', num2str(Times(i)), 'time' num2str(j) '.mat'];
        load(filename);
        if result2.info ~= 0
            disp(filename)
        end
        
        obj(i,j) = result2.obj/result2.params.N;
        u0(:,i,j) = result2.X(1:result2.params.ncontrols);
        Ks(:,i,j) = result2.X(end-1:end);
    end
end

meanu0 = mean(u0,3);
stdu0 = std(u0,[],3);
meanKs = mean(Ks,3);
stdKs = std(Ks,[],3);

figure
errorbar(Ns,meanu0(1,:), stdu0(1,:))
hold on;
errorbar(Ns,meanu0(2,:), stdu0(2,:), 'g')

figure;
errorbar(Ns,meanKs(1,:), stdKs(1,:))
hold on;
errorbar(Ns,meanKs(2,:), stdKs(2,:), 'g')
legend('Position', 'Derivative')

figure
errorbar(Ns,mean(obj,2), std(obj,[],2))