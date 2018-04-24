clear all
close all
clc

% load('1103Results_withcocon_big.mat'); stdevs = stdevs(1:end-3);
load('3103Results_better.mat');%('3003Results_Ubound.mat');%load('2803Results_Ton.mat');
Ws = Ws(1:end-1);
if ~isreal(lengthi)
    lengthi = 1;
end
params = result(1).result(1).result(1).result.params;

%Look at inputs
u = zeros(length(stdevs),params.ncontrols,length(Ws),lengthi);
Ks = zeros(length(stdevs),params.Ks,length(Ws),lengthi);
info = zeros(length(stdevs),length(Ws),lengthi);
obj = zeros(length(stdevs),length(Ws),lengthi);
cend = zeros(length(stdevs),length(Ws),lengthi);

for i = 1:lengthi
    for j = 1:length(stdevs)
        for k = 1:length(Ws)
            c = confun(result(j).result(k).result(i).result.X,result(j).result(k).result(i).result.params);
            cend(j,k,i) = c(end);
            info(j,k,i) = result(j).result(k).result(i).result.info;
            u(j,:,k,i) = result(j).result(k).result(i).result.X(1:2);
            obj(j,k,i) = result(j).result(k).result(i).result.obj;
            Ks(j,:,k,i) = result(j).result(k).result(i).result.X(end-3:end);
        end
    end
end

%Find minimum objective for each stdev and weight
minobj =zeros(size(info,1),size(info,2));
for i = 1:length(stdevs)
    for j = 1:length(Ws)
        minobj(i,j) = min(obj(i,j,:));
    end
end

meanobj =zeros(size(info,1),size(info,2));
stdobj =zeros(size(info,1),size(info,2));

meanu =zeros(size(u,1),size(u,2),size(u,3));
stdu =zeros(size(u,1),size(u,2),size(u,3));

meanK =zeros(size(Ks,1),size(Ks,2),size(Ks,3));
stdK =zeros(size(Ks,1),size(Ks,2),size(Ks,3));

for i = 1:size(info,1)
    for j = 1:size(info,2)
        objuse = [];
        for k = 1:size(info,3)
            if and(info(i,j,k) == 0, obj(i,j,k)<1.1*minobj(i,j))
                objuse = [objuse; obj(i,j,k)];
            end
        end
        meanobj(i,j) = mean(objuse);
        stdobj(i,j) = std(objuse);
    end
end

for i = 1:size(u,1)
    for j = 1:size(u,2)
        for h = 1:size(u,3)
            u_use = [];
            for k = 1:size(u,4)
                if and(info(i,j,k) == 0, obj(i,j,k)<1.1*minobj(i,j))
                    u_use = [u_use; u(i,j,h,k)];
                end
            end
            meanu(i,j,h) = mean(u_use);
            stdu(i,j,h) = std(u_use);
        end
    end
end
                
for i = 1:size(Ks,1)
    for j = 1:size(Ks,2)
        for h = 1:size(Ks,3)
            Ks_use = [];
            for k = 1:size(Ks,4)
                if and(info(i,j,k) == 0, obj(i,j,k)<1.1*minobj(i,j))
                    Ks_use = [Ks_use; Ks(i,j,h,k)];
                end
            end
            meanK(i,j,h) = mean(Ks_use);
            stdK(i,j,h) = std(Ks_use);
        end
    end
end


% figure; hold on;
% errorbar(stdevs, meanu(:,1,end),stdu(:,1,end),'Linewidth', 2.5)
% errorbar(stdevs, meanu(:,2,end),stdu(:,2,end),'Linewidth', 2.5)
% xlabel('Standard Deviation')
% ylabel('Co-contraction Input [-]')
% legend('Muscle 1', 'Muscle 2')
% set(gca,'Fontsize',16)
% 
% figure; hold on;
% errorbar(log(Ws)/log(10), meanu(end,1,:),stdu(end,1,:),'Linewidth', 2.5)
% errorbar(log(Ws)/log(10), meanu(end,2,:),stdu(end,2,:),'Linewidth', 2.5)
% xlabel('Weigth (logarithmic scale)')
% ylabel('Co-contraction Input [-]')
% legend('Muscle 1', 'Muscle 2')
% set(gca,'Fontsize',16)

dplotu(:,:) = meanu(:,1,:);
figure;
surf(log(Ws)/log(10),stdevs,dplotu)
xlabel('Weight')
ylabel('Standard Deviation')
zlabel('Open-loop Input u_0')
set(gca,'fontsize',14)

dplotu(:,:) = meanu(:,2,:);
figure;
surf(log(Ws)/log(10),stdevs,dplotu)
xlabel('Weight')
ylabel('Standard Deviation')
zlabel('Open-loop Input u_0')
set(gca,'fontsize',14)

figure; hold on;
errorbar(stdevs, meanK(:,1,end),stdK(:,1,end),'Linewidth', 2.5)
errorbar(stdevs, meanK(:,2,end),stdK(:,2,end),'Linewidth', 2.5)
xlabel('Standard Deviation')
ylabel('Position Feedback Gain')
legend('Muscle 1', 'Muscle 2')
set(gca,'Fontsize',16)

figure; hold on;
errorbar(log(Ws)/log(10), meanK(end,1,:),stdK(end,1,:),'Linewidth', 2.5)
errorbar(log(Ws)/log(10), meanK(end,2,:),stdK(end,2,:),'Linewidth', 2.5)
xlabel('Weigth (logarithmic scale)')
ylabel('Position Feedback Gain')
legend('Muscle 1', 'Muscle 2')
set(gca,'Fontsize',16)

figure; hold on;
errorbar(stdevs, meanK(:,3,end),stdK(:,3,end),'Linewidth', 2.5)
errorbar(stdevs, meanK(:,4,end),stdK(:,4,end),'Linewidth', 2.5)
xlabel('Standard Deviation')
ylabel('Derivative Feedback Gain')
legend('Muscle 1', 'Muscle 2')
set(gca,'Fontsize',16)

figure; hold on;
errorbar(log(Ws)/log(10), meanK(end,3,:),stdK(end,3,:),'Linewidth', 2.5)
errorbar(log(Ws)/log(10), meanK(end,4,:),stdK(end,4,:),'Linewidth', 2.5)
xlabel('Weigth (logarithmic scale)')
ylabel('Derivative Feedback Gain')
legend('Muscle 1', 'Muscle 2')
set(gca,'Fontsize',16)

%Look at RMS of targetangle
X = zeros(params.N,length(stdevs),length(Ws),1);
for i = 1:lengthi
    for j = 1:length(stdevs)
        for k = 1:length(Ws)
            Xuse = result(j).result(k).result(i).result.X;
            X1 = Xuse(params.ncontrols+1:end-params.Ks);
            x1 = reshape(X1, params.nstates, params.N, params.Nsamples);
            X(:,j,k,i) = x1(1,:);
        end
    end
end

Xerror(:,:,:) = rms((X-pi/2)/pi*180);
if lengthi == 1
    meanXerror(:,:) = Xerror(1,:,:);
else
    meanXerror(:,:) = mean(Xerror,3);
end

figure;
surf(log(Ws)/log(10),stdevs,meanXerror)
xlabel('Weigth (logarithmic scale)')
ylabel('Standard Deviation')
zlabel('RMS error from upright position [deg]')
