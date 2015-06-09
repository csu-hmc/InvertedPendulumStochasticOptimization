clear all
close all
clc

xx = [1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 110 130 150 170 190 200 300 400 500 600 700 800];

RMS = zeros(size(xx));
obj = zeros(size(xx));
K = zeros(2,length(xx));
meanx = zeros(length(xx),60);
for i = 1:length(xx)
    loadname = ['Solution_' num2str(xx(i)) 'SUs.mat'];
    load(loadname)
    obj(i) = result.obj;
    K(:,i) = result.X(end-1:end);
    [uperSU, uall,meanxx] = plotresult(result, 0);
    meanx(i,:) = meanxx(1,:);
    rmsu = rms(uall);
    RMS(i) = mean(rmsu(:));
end

figure
plot(xx,K')
figure
plot(xx,obj)
figure
plot(xx,RMS)