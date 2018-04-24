clear all
close all
clc

load('1103Results_withcocon_big.mat'); stdevs = stdevs(1:end-3);
params = result6(1).result(1).params;

%Look at inputs
u = zeros(length(stdevs),params.ncontrols,length(Ws),lengthi);
Ks = zeros(length(stdevs),params.Ks,length(Ws),lengthi);
info = zeros(length(stdevs),length(Ws),lengthi);
obj = zeros(length(stdevs),length(Ws),lengthi);
cend = zeros(length(stdevs),length(Ws),lengthi);

for i = 1:lengthi
%     for j = 2:length(stdevs)
        for k = 1:length(Ws)
            c = confun(result1(k).X,result1(k).params);
            cend(1,k,i) = c(end);
            info(1,k,i) = result1(k).info;
            u(1,:,k,i) = result1(k).X(1:2);
            obj(1,k,i) = result1(k).obj;
            Ks(1,:,k,i) = result1(k).X(end-3:end);
            
            c = confun(result2(k).result(i).X,result2(k).result(i).params);
            cend(1,k,i) = c(end);
            info(2,k,i) = result2(k).result(i).info;
            obj(2,k,i) = result2(k).result(i).obj;
            u(2,:,k,i) = result2(k).result(i).X(1:2);
            Ks(2,:,k,i) = result2(k).result(i).X(end-3:end);

            c = confun(result3(k).result(i).X,result3(k).result(i).params);
            cend(3,k,i) = c(end);
            info(3,k,i) = result3(k).result(i).info;
            obj(3,k,i) = result3(k).result(i).obj;
            u(3,:,k,i) = result3(k).result(i).X(1:2);
            Ks(3,:,k,i) = result3(k).result(i).X(end-3:end);

            c = confun(result4(k).result(i).X,result4(k).result(i).params);
            cend(4,k,i) = c(end);
            info(4,k,i) = result4(k).result(i).info;
            obj(4,k,i) = result4(k).result(i).obj;
            u(4,:,k,i) = result4(k).result(i).X(1:2);
            Ks(4,:,k,i) = result4(k).result(i).X(end-3:end);

            c = confun(result5(k).result(i).X,result5(k).result(i).params);
            cend(5,k,i) = c(end);
            info(5,k,i) = result5(k).result(i).info;
            obj(5,k,i) = result5(k).result(i).obj;
            u(5,:,k,i) = result5(k).result(i).X(1:2);
            Ks(5,:,k,i) = result5(k).result(i).X(end-3:end);

            c = confun(result6(k).result(i).X,result6(k).result(i).params);
            cend(6,k,i) = c(end);
            info(6,k,i) = result6(k).result(i).info;
            obj(6,k,i) = result6(k).result(i).obj;
            u(6,:,k,i) = result6(k).result(i).X(1:2);
            Ks(6,:,k,i) = result6(k).result(i).X(end-3:end);
 
%             c = confun(result7(k).result(i).X,result7(k).result(i).params);
%             cend(7,k,i) = c(end);
%             info(7,k,i) = result7(k).result(i).info;
%             obj(7,k,i) = result7(k).result(i).obj;
%             u(7,:,k,i) = result7(k).result(i).X(1:2);
%             Ks(7,:,k,i) = result7(k).result(i).X(end-3:end);
%             
%             c = confun(result8(k).result(i).X,result8(k).result(i).params);
%             cend(8,k,i) = c(end);
%             info(8,k,i) = result8(k).result(i).info;
%             obj(8,k,i) = result8(k).result(i).obj;
%             u(8,:,k,i) = result8(k).result(i).X(1:2);
%             Ks(8,:,k,i) = result8(k).result(i).X(end-3:end);
            
%             c = confun(result9(k).result(i).X,result9(k).result(i).params);
%             cend(9,k,i) = c(end);
%             info(9,k,i) = result9(k).result(i).info;
%             obj(9,k,i) = result9(k).result(i).obj;
%             u(9,:,k,i) = result9(k).result(i).X(1:2);
%             Ks(9,:,k,i) = result9(k).result(i).X(end-3:end);
        end
%     end
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
            if info(i,j,k) == 0
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
                if info(i,j,k) == 0
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
                if info(i,j,k) == 0
                    Ks_use = [Ks_use; Ks(i,j,h,k)];
                end
            end
            meanK(i,j,h) = mean(Ks_use);
            stdK(i,j,h) = std(Ks_use);
        end
    end
end


figure; hold on;
errorbar(stdevs, meanu(:,1,end),stdu(:,1,end),'Linewidth', 2.5)
errorbar(stdevs, meanu(:,2,end),stdu(:,2,end),'Linewidth', 2.5)
xlabel('Standard Deviation')
ylabel('Co-contraction Input [-]')
legend('Muscle 1', 'Muscle 2')
set(gca,'Fontsize',16)

figure; hold on;
errorbar(log(Ws)/log(10), meanu(end,1,:),stdu(end,1,:),'Linewidth', 2.5)
errorbar(log(Ws)/log(10), meanu(end,2,:),stdu(end,2,:),'Linewidth', 2.5)
xlabel('Weigth (logarithmic scale)')
ylabel('Co-contraction Input [-]')
legend('Muscle 1', 'Muscle 2')
set(gca,'Fontsize',16)

dplotu(:,:) = meanu(:,1,:);
figure;
surf(log(Ws)/log(10),stdevs,dplotu)
xlabel('Weight')
ylabel('Standard Deviation')
zlabel('Co-contraction Input u_0')
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
X = zeros(params.N,length(stdevs),length(Ws),lengthi);
for i = 1:lengthi
    for k = 1:length(Ws)
        Xuse = result1(k).X;
        X1 = Xuse(params.ncontrols+1:end-params.Ks);
        x1 = reshape(X1, params.nstates, params.N, params.Nsamples);
        X(:,1,k,i) = x1(1,:);
        Xuse = result2(k).result(i).X;
        X1 = Xuse(params.ncontrols+1:end-params.Ks);
        x1 = reshape(X1, params.nstates, params.N, params.Nsamples);
        X(:,2,k,i) = x1(1,:);
        Xuse = result3(k).result(i).X;
        X1 = Xuse(params.ncontrols+1:end-params.Ks);
        x1 = reshape(X1, params.nstates, params.N, params.Nsamples);
        X(:,3,k,i) = x1(1,:);
        Xuse = result4(k).result(i).X;
        X1 = Xuse(params.ncontrols+1:end-params.Ks);
        x1 = reshape(X1, params.nstates, params.N, params.Nsamples);
        X(:,4,k,i) = x1(1,:);
        Xuse = result5(k).result(i).X;
        X1 = Xuse(params.ncontrols+1:end-params.Ks);
        x1 = reshape(X1, params.nstates, params.N, params.Nsamples);
        X(:,5,k,i) = x1(1,:);
        Xuse = result6(k).result(i).X;
        X1 = Xuse(params.ncontrols+1:end-params.Ks);
        x1 = reshape(X1, params.nstates, params.N, params.Nsamples);
        X(:,6,k,i) = x1(1,:);
        Xuse = result7(k).result(i).X;
        X1 = Xuse(params.ncontrols+1:end-params.Ks);
        x1 = reshape(X1, params.nstates, params.N, params.Nsamples);
        X(:,7,k,i) = x1(1,:);
        Xuse = result8(k).result(i).X;
        X1 = Xuse(params.ncontrols+1:end-params.Ks);
        x1 = reshape(X1, params.nstates, params.N, params.Nsamples);
        X(:,8,k,i) = x1(1,:);
    end
end

Xerror(:,:,:) = rms((X-pi/2)/pi*180);
meanXerror(:,:) = mean(Xerror,3);

figure;
surf(log(Ws)/log(10),stdevs,meanXerror)
xlabel('Weigth (logarithmic scale)')
ylabel('Standard Deviation')
zlabel('RMS error from upright position [deg]')
