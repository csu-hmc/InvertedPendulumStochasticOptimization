function [Fsee1,Fsee2, Fsee3] = testMuscDyns(params)

    nmus = params.nmus;
    ndof = params.ndof;
    x = [pi/2;0;1;1;1;1];
    x_org = x;
    xdot = zeros(size(x));
    xdot_org = xdot;

    %Test lce
    lcearr = 0.5:0.05:1.5;
    Fsee1 = zeros(nmus,length(lcearr));
    for i = 1:nmus
        for j = 1:length(lcearr)
            x(ndof*2+nmus+i) = lcearr(j);
            Fsee = getMusDyns(x,xdot,[0;0],params);
            Fsee1(i,j) = Fsee(i);
        end
        x = x_org;
    end

    %test act
    actarr = 0:0.05:1;
    Fsee2 = zeros(nmus,length(actarr));
    for i = 1:nmus
        for j = 1:length(lcearr)
            x(ndof*2+i) = actarr(j);
            Fsee = getMusDyns(x,xdot,[0;0],params);
            Fsee2(i,j) = Fsee(i);
        end
        x = x_org;
    end
    
    %test lcedot
    lcedotarr = -10:0.5:10;
    Fsee3 = zeros(nmus,length(lcedotarr));
    for i = 1:nmus
        for j = 1:length(lcedotarr)
            xdot(ndof*2+nmus+i) = lcedotarr(j);
            Fsee = getMusDyns(x,xdot,[0;0],params);
            Fsee3(i,j) = Fsee(i);
        end
        xdot = xdot_org;
    end
    
    figure
    subplot(3,1,1)
    plot(lcearr,Fsee1)
    xlabel('Lce')
    ylabel('Flce')
    
    subplot(3,1,2)
    plot(actarr,Fsee2)
    xlabel('Act')
    ylabel('Flce')

    subplot(3,1,3)
    plot(lcedotarr,Fsee3)
    xlabel('lcedot')
    ylabel('Flce')

end