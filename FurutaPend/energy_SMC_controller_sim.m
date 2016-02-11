% MCE 647
% Final Project
% Holly Warner

function out = energy_SMC_controller_sim(z)

% parse inputs
q1 = z(1);
q2 = z(2);
q1dot = z(3);
q2dot = z(4);

% Swing up energy controller tuning
E0 = 0;
k_E = 40.75;

% Stabilization sliding mode controller tuning
C = [1; 3000; 1; 750];
a = .2;
k_SMC = 2;

% Switching parameters
E_threshold = .1;
q2_threshold = .2;

% define physical system parameters
m2 = 0.1022 + .0448; %kg
I2z = 1/12*m2*((3/4*2.54/100)^2+(10*2.54/100)^2); %kg-m^2
lc2 = 5*2.54/100; %m

g = 9.81; %m/s^2

% calculate pendulum total energy
w0 = sqrt(m2*g*lc2/I2z);
E = m2*g*lc2*(1/2*(q2dot/w0)^2+cos(q2)-1);

% check if threshold requirements are met
% true: sliding mode controller
% false: swing-up controller
if (E > -E_threshold && E < E_threshold) && ...
   ((q2 > -q2_threshold && q2 < q2_threshold) || ...
   (q2 > (2*pi-q2_threshold) && q2 < (2*pi+q2_threshold)))
    
    % shift incremental encoder data
    if (q2 > (2*pi-q2_threshold)) && (q2 < (2*pi+q2_threshold))
        q2 = q2-2*pi;
    else
        q2 = q2;
    end
    
    % linear state space model
    A = [0   0                   1                  0;
        0   0                   0                  1;
        0   0.266410130501494  -3.440732130762066  0;
        0  58.250770915294339  -5.153850593745645  0];
    
    B = [0; 0; 4.529374727031972; 0];
    
    % error vector
    q_qdot_tilde = [q1;q2;q1dot;q2dot]-[0;0;0;0];
    
    % sliding surface
    s = C'*q_qdot_tilde;
    
    % sliding mode control law
    u = -inv(C'*B)*(C'*A*[q1;q2;q1dot;q2dot] + k_SMC*abs(s)^a*sign(s));
    
else
    
    % swing-up energy control law
    u = -(k_E*(E-E0))*sign(q2dot*cos(q2));
    
    % default s to distinguish between controllers
    s=-1;
    
end

% output final control signal in Amps, pendulum total energy, and s
out = [u/4.664;E;s];

end