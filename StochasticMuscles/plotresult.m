function [x1,u0, uall] = plotresult(result, plotyes)

if nargin == 1
    plotyes = 1;
end
params = result.params;
if ~isfield(params, 'Ks')
    params.Ks = 2;
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
        Fsee(:,i,j) = getMusDyns_asat(x1(:,i,j),zeros(6,1),uall(:,i,j),params);
        Torque(i,j) = params.muscleparam.d*Fsee(:,i,j);
    end
end

Fseem = mean(Fsee,3);
Torquem = mean(Torque,2);

if plotyes == 1
    T = result.params.T;
    h = result.params.h;
    figure
    hold on
    for i = 1:size(x1,3)
        plot([0:h:T],x1(1,:,i)/pi*180,'color',0.2*[1 1 1])
    end
    plot([0:h:T],meanx(1,:)/pi*180, 'b', 'LineWidth', 1.5);
    plot([0 T], [pi/2 pi/2]/pi*180, 'r')
    xlabel('Time [s]')
    ylabel('Pendulum Angle [deg]')
    
    figure
%     plot([0:h:T],umean)
%     hold on
%     plot([0:h:T],uall(:,:,1), '--')
    plot(transpose(uall(:,:)), '--')
    xlabel('Time [s]')
    ylabel('Input signal')
    
    figure
    plot([0:h:T], Torquem)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    
    figure
    subplot(3,1,1)
%     plot([0:h:T],Fseem')
    plot(transpose(Fsee(:,:)))
    ylabel('Muscle force Fsee [N]')
    
    subplot(3,1,2)
    plot([0:h:T],meanx(3:4,:))
    plot(transpose(x1(3:4,:)))
    ylabel('Activation State')
    
    subplot(3,1,3)
%     plot([0:h:T],meanx(5:6,:))
    plot(transpose(x1(5:6,:)))
    xlabel('Time [s]')
    ylabel('CE length [-]')
    
end