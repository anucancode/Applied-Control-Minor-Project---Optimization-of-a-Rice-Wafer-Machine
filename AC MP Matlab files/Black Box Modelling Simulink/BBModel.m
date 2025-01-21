clear all
close all
clc

%% Import Data
load('BBModel.mat');
load('long measured data interpolated.mat')
plot(tmin,UpperCookingPlate,'linewidth',3);
ylabel('Temperature [T]');
xlabel('time [min]');
title('Cleaned Data');

P_bc_cleaned=P_bc(:,2);
T_bhp_measured_interp_cleaned=T_bhp_measured_interp(:,2);

%%
systemIdentification
