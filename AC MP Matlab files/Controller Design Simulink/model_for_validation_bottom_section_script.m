%% Prepare workspace
clear all
close all
clc

%% Parameters
diameter_cilinder = 0.092; % diameter of cilinder [m]
area_cilinder = pi * (diameter_cilinder/2)^2; % area of cilinder [m2]

% top metal
density_m = 7700; % density [kg/m^3]
c_tm = 460; % specific heat capacity [J/(kgC)]
k_tm = 22; % thermal conductivity [W/(mC)]
d_tm = 0.035; % thickness [m]
m_tm = area_cilinder*d_tm*density_m; % mass [kg]
A_tm = area_cilinder; % cross area [m^2]
A_tm_air = pi*diameter_cilinder*d_tm; % outside area [m^2]

% bottom metal
m_bm = m_tm; % mass [kg]
c_bm = c_tm; % specific heat capacity [J/(kgC)]
k_bm = k_tm; % thermal conductivity [W/(mC)]
d_bm = d_tm; % thickness [m]
A_bm = area_cilinder; % cross area [m^2]
A_bm_air = A_tm_air; % outside area [m^2]

% coil
% +- 10% resistance
% +- 2.5% heated length
% coil parameters taken from datasheet, crni steel seems to have similar
% properties as the other metal (sources)
density_c = density_m; % density [kg/m^3]
length_c = 0.9; % length [m]
w_c = 0.0042; % coil width [m]
h_c = 0.0022; % coil height [m]
c_tc = c_tm;
d_tc = h_c;
m_tc = w_c*h_c*length_c*density_c; % mass coil [kg]
d_bc = d_tc;
c_bc = c_tc;
m_bc = m_tc;

% top isolation
density_i = 2200; % density [kg/m^3]
c_ti = 880; % specific heat capacity [J/(kgC)]
k_ti = 0.3; % thermal conductivity [W/(mC)]
d_ti = 0.02; % thickness [m]
m_ti = area_cilinder*d_ti*density_i; % mass [kg]
A_ti = area_cilinder; % cross area [m^2]
A_ti_air = pi*diameter_cilinder*d_ti; % outside area [m^2]

% bottom isolation
m_bi = m_ti; % mass [kg]
c_bi = c_ti; % specific heat capacity [J/(kgC)]
k_bi = k_ti; % thermal conductivity [W/(mC)]
d_bi = d_ti; % thickness [m]
A_bi = area_cilinder; % cross area [m^2]
A_bi_air = A_ti_air; % outside area [m^2]

% top heating plate
c_thp = c_tm; % specific heat capacity [J/(kgC)]
k_thp = k_tm; % thermal conductivity [W/(mC)]
d_thp = 0.01; % thickness [m]
m_thp = area_cilinder*d_thp*density_m; % mass [kg]
A_thp = area_cilinder; % cross area [m^2]
A_thp_air = pi*diameter_cilinder*d_thp; % outside area [m^2]

% bottom heating plate
m_bhp = m_thp; % mass [kg]
c_bhp = c_tm; % specific heat capacity [J/(kgC)]
k_bhp = k_tm; % thermal conductivity [W/(mC)]
d_bhp = d_thp; % thickness [m]
A_bhp = area_cilinder; % cross area [m^2]
A_bhp_air = A_thp_air; % outside area [m^2]

% ring
density_r = 7800; % density [kg/m^3]
c_r = 470; % specific heat capacity [J/(kgC)]
k_r = 22.8; % thermal conductivity [W/(mC)]
d_r = 0.0075; % ring thickness [m]
h_r = 0.0245; % ring height [m]
diam_ring_outer = diameter_cilinder+2*d_r; % outer ring diameter [m]
area_ring_outer = pi * (diam_ring_outer/2)^2;
h_r_tm = h_r-d_thp; % height of ring which overlaps with the metal [m]
m_r = pi*(diam_ring_outer/2)^2*h_r*density_r-area_cilinder*h_r*density_r; % mass ring [kg]
A_r_air = pi*diam_ring_outer*h_r + 2 * area_ring_outer - area_cilinder; % area of ring exposed to air [m^2]
A_r_tm = pi*diameter_cilinder*h_r_tm; % area of ring overlapping with metal in top position [m^2]

% air
h_c = 26.5; % convective heat transfer coefficient [W/(m^2C)] (see source pdf)
T_a = 18;

%% Set up controller and disturbance
% controller settings
P = 100; 
I = 0;
D = 0;

T_setpoint = 285; % temperature setpoint

% baking for 3 seconds every 12 seconds, total 700W extra compared to
% steady state, so assume 2800W extracted during 3 secs to account for the
% extra inputted wattage = 1400W per plate for 3 secs
% to check wattage used
P_rice = 1400;
T_init = T_a; % initial conditions (for all temperatures!)


%% Simulate controller
out = sim('model_for_validation_bottom_section.slx', 'StartTime','0','StopTime','7200','FixedStep','0.01');

%out = sim('model_for_validation_bottom_section_mpc.slx', 'StartTime','0','StopTime','7200','FixedStep','0.01');
%% Plot results
close all
figure;
plot(out.T_bhp);
hold on
plot(out.T_bc);
plot(out.T_bm);
plot(out.T_bi);
plot(out.T_bhp.Time, abs(out.T_bhp.Data - T_setpoint));
title("Bottom Section Simulated Temperatures")
ylabel("Temperature [C]")
yyaxis right
plot(out.P_bc);
%plot(out.P_rice);
ylabel("Coil Power [W]")
legend("T thp", "T tc", "T tm", "T ti", "error", "P bc", 'location', 'e');