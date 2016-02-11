function pendanim(x,params, moviefile)
% pendanim.m: make a movie of the pendulum swing up		

	% Vector of time instants for nodes
	T = params.T; 
	
	% initialize movie file
    if nargin == 2
        moviefile = 'pendsu.avi';
    end
    
	avi = VideoWriter(moviefile);% 'fps', 30, 'compression', 'xvid','quality',95);
    avi.FrameRate = 15;
    
	% initialize figure window
	close all
	figure(1);
	clf;
	set(gcf,'Position',[5 5 550 550]);
	set(gcf, 'color', 'white');
	
% 	% create mark to indicate zero?
% 	np = 15000;
% 	xg = rand(np,1);
% 	yg = rand(np,1);
% 	xg = -1 + 3*[xg ; 2-xg];
% 	yg = -0.15*[yg ; yg];
	
    l = params.l; %length of the pendulum
    
	% make the movie
	R = [1:6 4];			% right stick points
	L = [2 7:10 8];			% left stick points
	nframes = params.NperSU;
    %Always start at zero position
    xx(1) = 0;
    yy(1) = 0;
    open(avi);
    for i=0:nframes-1
%         plot(xg,yg,'.','Color',[0.7 0.7 0.7],'MarkerSize',4);
%         hold on

        % Find positions
        xx(2) = l*cos(x(1,i+1));
        yy(2) = l*sin(x(1,i+1));
        plot(xx,yy,'b','LineWidth',2);
        axis([-1.5 1.5 -1.5 1.5]);
        axis('off');
        if (i==0)
            F = getframe(gca);
            frame = [1 1 size(F.cdata,2) size(F.cdata,1)];
        else
            F = getframe(gca,frame);
        end
        writeVideo(avi,F);
        cla;
	end
	close(avi);
	hold off;
	close all
end