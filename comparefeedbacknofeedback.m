function comparefeedbacknofeedback(fresult, nfresult)

%First do things for with feedback
params = fresult.params;
X = fresult.X;
x = reshape(X(1:end-2),params.nvarpernode, params.N);
K = X(end-1);
Kd = X(end);

xperSU = zeros(params.nstates, params.NperSU, params.NSU);
uperSU = zeros(params.ncontrols,params.NperSU, params.NSU);
for i = 1:params.NSU
    xperSU(:,:,i) = x(1:params.nstates,(i-1)*params.NperSU+(1:params.NperSU));
    uperSU(:,:,i) = x(params.nstates+1:end,(i-1)*params.NperSU+(1:params.NperSU));
    uperSU(:,:,i) = uperSU(:,:,i) + K*xperSU(1,:,i) + Kd*xperSU(2,:,i);
end

meanxf = mean(xperSU,3);
stdxf = std(xperSU,[],3);
meanuf = mean(uperSU,3);


%Now do things without feedback
params = nfresult.params;
X = nfresult.X;
x = reshape(X(1:end-2),params.nvarpernode, params.N);

xperSU = zeros(params.nstates, params.NperSU, params.NSU);
uperSU = zeros(params.ncontrols,params.NperSU, params.NSU);
for i = 1:params.NSU
    xperSU(:,:,i) = x(1:params.nstates,(i-1)*params.NperSU+(1:params.NperSU));
    uperSU(:,:,i) = x(params.nstates+1:end,(i-1)*params.NperSU+(1:params.NperSU));
end

meanxnf = mean(xperSU,3);
stdxnf = std(xperSU,[],3);
meanunf = mean(uperSU,3);

%Plotting is awesome of course
T = fresult.params.T;
h = T/(fresult.params.NperSU-1);
% xstds = [0:h:T T:-h:0];
% plotstds = meanx-stdx;
% plotstds(:,end+1:2*end) = meanx(:,end:-1:1)+stdx(:,end:-1:1);
figure
hold on
plot([0:h:T],meanxf(1,:), 'b');
plot([0:h:T],meanxnf(1,:), 'k');
plot([0 T], [pi/2 pi/2], 'r')
figure
plot([0:h:T],meanuf(1,:), 'b')
hold on
plot([0:h:T],meanunf(1,:), 'k')