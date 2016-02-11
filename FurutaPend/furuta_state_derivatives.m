% MCE 647
% Final Project
% Holly Warner

function zdot = furuta_state_derivatives(z,i_a)

% parse inputs
q1 = z(1);
q2 = z(2);
q1dot = z(3);
q2dot = z(4);

% define physical parameters
m1 = 1.2706; %kg
m2 = 0.1022+.0156; %kg
I1y = 1/12*m1*((1/2*2.54/100)^2+(10*2.54/100)^2); %kg-m^2
I2z = 1/12*m2*((3/4*2.54/100)^2+(10*2.54/100)^2); %kg-m^2
l1 = 10*2.54/100; %m
lc2 = 5*2.54/100; %m

g = 9.81; %m/s^2
 
% inertia matrix
D = [I1y + (l1^2*m1)/4 + l1^2*m2 + lc2^2*m2*sin(q2)^2 + 1,...
    -l1*lc2*m2*cos(q2);...
    -l1*lc2*m2*cos(q2),...
    m2*lc2^2 + I2z];
 
% coriolis matrix
C = [(lc2^2*m2*q2dot*sin(2*q2))/2, lc2*m2*sin(q2)*...
    (l1*q2dot + lc2*q1dot*cos(q2));...
    -(lc2^2*m2*q1dot*sin(2*q2))/2,...
    0];
 
% losses matrix
R = [3.543 0; 0 0];

% gravity vector
gq = [0; -g*lc2*m2*sin(q2)];

% system input matrix
u = [4.664*i_a; 0];

% state derivatives
zdot = [q1dot;q2dot;inv(D)*(u-C*[q1dot;q2dot]-R*[q1dot;q2dot]-gq)];

end