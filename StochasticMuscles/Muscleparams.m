function muscleparam = Muscleparams()

muscleparam.A = 0.25;
muscleparam.vmax = 10;                      %-
muscleparam.w = 0.56;
muscleparam.Ts = [0.0120 0.0476];           %s[0.120 0.476]; %
muscleparam.d = [0.02 -0.02];% [0.024 -0.02];              %m
muscleparam.fmax = [1100 1100]; % [1039 1144];             %N
muscleparam.gmax = 1.5;
muscleparam.names = ['brachialis', 'triceps'];    
muscleparam.lceopt = [0.07 0.07]; % [0.0738 0.0691];       %m
muscleparam.lsees = [0.05 0.05]; %[0.0392 0.1706];        %m
muscleparam.lpees = [1.2 1.2];              %-
muscleparam.l0 = muscleparam.lsees+muscleparam.lceopt;% [0.13774 0.21641];        %m
muscleparam.kpee = [1 1]; 
muscleparam.ksee = 1./(0.05^2*muscleparam.lsees.^2); %from gait2dc




