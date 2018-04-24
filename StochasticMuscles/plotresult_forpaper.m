function [x1,u0, uall] = plotresult_forpaper(result, result1)

if nargin == 0
    load('051216PeriodicResults.mat');
    result = result6(end).result(end);
    result1 = result1(end);
end

% result
params = result.params;
if ~isfield(params, 'Ks')
    paramsks = 2;
end
X = result.X;
X1 = X(params.ncontrols+1:end-params.Ks);
x1 = reshape(X1, params.nstates, params.N, params.Nsamples);
u0 = X(1:params.ncontrols);

uall = zeros(params.nmus,params.N,params.Nsamples);
% Get complete input
for j = 1:params.Nsamples
    for i = 1:params.N
        uall(:,i,j) = findTorque(u0,X(end-params.Ks+1:end),x1(1:params.ndof*2,i,j), params);
    end
end

meanx = mean(x1,3);
umean = mean(uall,3);

Fsee = zeros(params.nmus,params.N,params.Nsamples);
Torque = zeros(params.N, params.Nsamples);
for j = 1:params.Nsamples
    for i = 1:params.N
        Fsee(:,i,j) = getMusDyns(x1(:,i,j),zeros(6,1),uall(:,i,j),params);
        Torque(i,j) = params.muscleparam.d*Fsee(:,i,j);
    end
end

Fseem = mean(Fsee,3);
Torquem = mean(Torque,2);

% result1
params = result1.params;
if ~isfield(params, 'Ks')
    paramsks = 2;
end
X1 = result1.X;
X11 = X1(params.ncontrols+1:end-params.Ks);
x11 = reshape(X11, params.nstates, params.N, params.Nsamples);
u01 = X1(1:params.ncontrols);

uall1 = zeros(params.nmus,params.N,params.Nsamples);
% Get complete input
for j = 1:params.Nsamples
    for i = 1:params.N
        uall1(:,i,j) = findTorque(u01,X1(end-params.Ks+1:end),x11(1:params.ndof*2,i,j), params);
    end
end

meanx1 = mean(x11,3);

Fsee1 = zeros(params.nmus,params.N,params.Nsamples);
Torque1 = zeros(params.N, params.Nsamples);
for j = 1:params.Nsamples
    for i = 1:params.N
        Fsee1(:,i,j) = getMusDyns(x11(:,i,j),zeros(6,1),uall1(:,i,j),params);
        Torque1(i,j) = params.muscleparam.d*Fsee1(:,i,j);
    end
end

Torquem1 = mean(Torque1,2);

T = result1.params.T;
h = result1.params.h;
figure
subplot(6,2,[1,3,5])
hold on
plot([0:h:T],meanx1(1,:)/pi*180, 'k', 'LineWidth', 1.2);
plot([0:h:T],meanx(1,:)/pi*180, 'r', 'LineWidth', 1.2);
% plot([0 T], [pi/2 pi/2]/pi*180, 'b')
xlabel('Time [s]')
ylabel('Pendulum Angle [deg]')
legend('Deterministic', 'Stochastic')

subplot(6,2,[7,9,11])
hold on
plot([0:h:T], transpose(uall(:,:)), 'r', 'LineWidth', 1.2)
plot([0:h:T], transpose(uall1(:,:)), 'k','LineWidth', 1.2)
legend('Muscle 1', 'Muscle 2')
xlabel('Time [s]')
ylabel('Input signal')

subplot(6,2,[8,10,12])
hold on
plot([0:h:T], Torquem, 'r','LineWidth', 1.2)
plot([0:h:T], Torquem1, 'k', 'LineWidth', 1.2)
xlabel('Time [s]')
ylabel('Torque [Nm]')

subplot(6,2,2)
hold on
plot([0:h:T],transpose(Fsee1(:,:)), 'k', 'LineWidth', 1.2)
plot([0:h:T],transpose(Fsee(:,:)), 'r', 'LineWidth', 1.2)
ylabel({'Muscle force';'Fsee [N]'})

subplot(6,2,4)
hold on
plot([0:h:T],transpose(x11(3:4,:)), 'k', 'LineWidth', 1.2)
plot([0:h:T],transpose(x1(3:4,:)), 'r', 'LineWidth', 1.2)
ylabel({'Activation'; 'State'})

subplot(6,2,6)
hold on
plot([0:h:T],transpose(x11(5:6,:)), 'k', 'LineWidth', 1.2)
plot([0:h:T],transpose(x1(5:6,:)), 'r', 'LineWidth', 1.2)
xlabel('Time [s]')
ylabel({'CE length'; '[-]'})
end