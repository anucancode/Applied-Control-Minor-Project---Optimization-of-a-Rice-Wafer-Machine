% Preparing the workspace
close all
clear all
clc

% Define the thermal model parameters
R = 1;  % Thermal resistance
C = 1;  % Thermal capacitance

% Define the desired temperature set point
setPoint = 100;

% Define the PID controller gains 
Kp = 1;
Ki = 0.5;
Kd = 0.2;

% Create the PID controller object
pidController = pid(Kp, Ki, Kd);

% Set the controller sample time 
sampleTime = 0.1;  % To be adjusted
pidController.Ts = sampleTime;

% Define the output limits for the PID controller
minOutput = -1;  % Minimum control signal value
maxOutput = 1;   % Maximum control signal value

% Simulate the rice popping process with the PID controller
tFinal = 10;  % Simulation time
t = 0:sampleTime:tFinal;
numSamples = numel(t);
temperature = zeros(numSamples, 1);
temperature(1) = 18;  % Initial temperature (T_a)

for i = 2:numSamples
    % Compute the control signal from the PID controller
    if i == 2
        controlSignal = pidController.Kp * (setPoint - temperature(i-1)) 
            + pidController.Ki * pidController.InputDelay * sum(setPoint - temperature(1:i-1)) 
            + pidController.Kd * (temperature(i-1) - temperature(1)) / sampleTime;
    else
        controlSignal = pidController.Kp * (setPoint - temperature(i-1)) 
            + pidController.Ki * pidController.InputDelay * sum(setPoint - temperature(1:i-1)) 
            + pidController.Kd * (temperature(i-1) - temperature(i-2)) / sampleTime;
    end
    
    % Apply output limits to the control signal
    controlSignal = max(minOutput, min(maxOutput, controlSignal));
    
    % Simulate the rice popping process dynamics
    temperatureChange = (controlSignal - R * (temperature(i-1) - setPoint)) / C;
    temperature(i) = temperature(i-1) + sampleTime * temperatureChange;
end

% Plot the temperature response
figure;
plot(t, temperature);
hold on;
plot([0 tFinal], [setPoint setPoint], 'r--');
xlabel('Time');
ylabel('Temperature');
legend('Temperature Response', 'Set Point');
title('PID Controller for Rice Popping Machine Thermal Model');






