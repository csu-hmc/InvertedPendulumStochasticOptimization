function isometrictest(params)
    disp(' ');
    disp('Calculating isometric strength curves...');
    clf;
    
    % hip moment-angle curves, at 30-deg intervals
    isometric_curves(1,[70:1:110],params); 
end

function isometric_curves(joint1, range1,params)
    angvel = 0;
    pascurves = [];
    poscurves = [];
    poslce = [];
    neglce = [];
    posF = [];
    negF = [];
    negcurves = [];
    for angle1 = range1+1e-6
        pascurves = [pascurves maxmoment(joint1, angle1, angvel, 0,params)];
        [mom, lce,Fsee] = maxmoment(joint1, angle1, angvel, 1,params);
        poscurves = [poscurves mom];
        poslce = [poslce lce'];
        posF = [posF Fsee];
        [mom, lce,Fsee] = maxmoment(joint1, angle1, angvel, -1,params);
        negcurves = [negcurves mom];
        neglce = [neglce lce'];
        negF = [negF Fsee];
    end

    % plot total moments on left side of figure
    subplot(3,1,3*joint1-2);
    plot(range1, poscurves);hold on;
    plot(range1, negcurves);
    labels;
    title(['Total moment']);

    % plot passive moments in middle column of figure
    subplot(3,1,3*joint1-1);
    plot(range1, pascurves);hold on;
    labels;
    title(['Passive moment']);

    % subtract passive moments and plot in rightmost column of figure
    subplot(3,1,3*joint1);
    plot(range1, poscurves-pascurves);hold on;
    plot(range1, negcurves-pascurves);
    labels;
    title(['Active = total -- passive']);
    
    figure
    subplot(2,1,1)
    plot(range1, poslce); hold on; plot(range1,neglce)
    xlabel('Angles [deg]')
    ylabel('Lce [m]')
    
    subplot(2,1,2)
    plot(range1, posF); hold on; plot(range1,negF)
    xlabel('Angles [deg]')
    ylabel('Fsee [N]')
    
%     legend(legends);
    %======================================================================
    %> @brief Adds labels to plot.
    %======================================================================
    function labels
        a = get(gca);
        axis([min(range1) max(range1) a.YLim]);
        plot([0 0],a.YLim,'k:');
        plot(a.XLim,[0 0],'k:');
        hold off;
        xlabel('angle (deg)');
        ylabel('moment (Nm)');
    end
end

function [mom,Lceans,Fsee] = maxmoment(joint, angles, angvel, sign,params)
    angles = angles*pi/180;		% convert to radians
    angvel = angvel*pi/180;

    % determine moment arms so we know which muscles to activate
    % here (2D model) we have constant moment arms.
    % we should in a later version ask the MEX function what the moment arms are at this posture
    momentarms = params.muscleparam.d';
    
    Act =  sign*momentarms(:,joint) > 0;	% vector that has a 1 for all muscles we want to activate

    % determine lengthening velocities of the muscle-tendon complexes, normalize to Lceopt
    Vmuscle = -(momentarms * angvel) ./ params.muscleparam.lceopt';

    % determine the Lce's at which there is contraction equilibrium (dF/dt = 0, or Lcedot = Vmuscle)
    x = [angles; angvel; Act; zeros(params.nmus,1)];
    xdot = [zeros(2*params.ndof,1) ; Vmuscle ; zeros(params.nmus,1)];	% we want these state derivatives
    u = zeros(params.nmus,1);				% no stim, we don't care about activation dynamics here

    % solve the equilibrium equation with Matlab fzero function, one muscle at a time
    for imus=1:params.nmus;
        Lce = 1.0;		% initial guess for this muscle's Lce			
        [Lceans(imus), Fval, Flag] = fzero(@contraction_equilibrium, Lce);				
    end

    if (flag < 0)
        fprintf('maxmoment: muscle contraction equilibrium not found within max number of iterations.\n');
        keyboard
    end

    % now determine the joint moments at this state of the system
    Fsee = getMusDyns(x,xdot,u,params);
    mom = params.muscleparam.d*Fsee;

    %======================================================================
    %> @brief Test dynamics.
    %>
    %> @param Lce
    %> @retval F
    %======================================================================
    function [F] = contraction_equilibrium(Lce)
        x(2*params.ndof+params.nmus+imus) = Lce;
        f = StocDyn(x, xdot, u, [0;0], params);
        F = f(2*params.ndof+params.nmus+imus);
    end

end