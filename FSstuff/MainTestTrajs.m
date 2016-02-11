clear all
close all
clc

% Do same thing using stochastic trajectory
load('Solutions\LinFeed\Solution_30SUs2.mat')
% load('TryForLQR.mat')
% load('WithManyK\WithManyKresult1.mat')

ratios = 1e6;%[ 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10]; %

xend = zeros(2,2,50,length(ratios));%);%
rmsu = zeros(2,50,length(ratios));%);%
rmsu0 = zeros(2,50,length(ratios));%);%

for j = 1:length(ratios)
    for i = 1:50
        s = rng;
        [xend(:,:,i,j), rmsu(:,i,j), rmsu0(:,i,j),K,K1] = getLQRresults(result(1),result(7),ratios(j),1);
    end
end

rmsx = rms(xend(:,1,:,:)-pi/2,3);
meantheta(:,:) = rmsx;
rmsdx = rms(xend(:,2,:,:),3);
meandtheta(:,:) = rmsdx;

meanu(:,:) = mean(rmsu,2);
meanu0(:,:) = mean(rmsu0,2);

figure
semilogx(ratios,meantheta')
xlabel('Q/R Ratio')
ylabel('Accuracy of Final State [rad]')

figure
semilogx(ratios,meandtheta')
xlabel('Q/R Ratio')
ylabel('Final Velocity [rad/s]')

figure
semilogx(ratios,meanu')
xlabel('Q/R Ratio')
ylabel('Root Mean Square Input [Nm]')

figure
semilogx(ratios,meanu'-meanu0')
xlabel('Q/R Ratio')
ylabel('Root Mean Square Feedback Control')