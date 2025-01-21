close all
clear all
clc
% Step 1: Obtain the transfer function of the PI controller
Kp = 1; % Proportional gain
Ti = 1; % Integral time constant
s = tf('s');
PI = Kp + (Kp/Ti)/s;

% Step 2: Determine the controller parameters (proportional gain and integral time constant)
Kp = 1.2; % Proportional gain
Ti = 2.5; % Integral time constant
PI = Kp + (Kp/Ti)/s;

% Step 3: Combine the controller transfer function with the process transfer function to form the open-loop transfer function
H = ((0.0030*s^3)+(4.8978e-5*s^2)+(2.5849e-7*s)+(4.3592e-10))/((s^4)+(0.0288*s^3)+(2.3163e-4*s^2)+(5.8211e-7*s)+(2.4497e-10));
G = PI * H;

% Step 4: Use the rlocus function to plot the root locus
rlocus(G);

% Step 5: Analyze the root locus plot to determine the desired controller parameters
% Determine the desired pole locations and corresponding gains from the root locus plot

% Step 6: Implement the PI controller in code
t = 0:0.01:10; % Time vector
r = ones(size(t)); % Reference input
[y, t] = lsim(feedback(G,1), r, t); % Simulate the closed-loop response

% Plot the response
plot(t, y);
xlabel('Time');
ylabel('Output');
title('Response of the Rice Popping Process with PI Controller (Root Locus)');
