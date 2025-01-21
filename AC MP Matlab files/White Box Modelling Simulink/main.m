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
m_r = area_ring_outer*h_r*density_r-area_cilinder*h_r*density_r; % mass ring [kg]
A_r_air = pi*diam_ring_outer*h_r + 2 * area_ring_outer - area_cilinder; % area of ring exposed to air [m^2]
A_r_tm = pi*diameter_cilinder*h_r_tm; % area of ring overlapping with metal in top position [m^2]

% air
% todo split into 26.5 for bot section, 16 for top section with ring
% try average too
h_c_top = 16.5; % convective heat transfer coefficient [W/(m^2C)] (see source pdf)
h_c_bot = 26.5;

%% validation setup for comparing to warmup in masters information 26-4
% power data taken from graph 1 final report jasper, from excel sheet
t = 0:60:900;
P_tc_org = [0 600 660 600 660 600 600 660 600 660 360 360 240 300 60 180];
P_tc = [[0:0.01:900]', interp1(t, P_tc_org', [0:0.01:900]')]; % top coil heat added [W]
P_bc_org = [0 660 600 660 600 660 660 600 420 300 240 300 180 240 180 60];
P_bc = [[0:0.01:900]', interp1(t, P_bc_org', [0:0.01:900]')]; % bottom coil head added [W]
P_rice = [0 0; 5400 0]; % rice heat consumed [W]
T_a = 18; % ambient temperature [C]
T_init = T_a; % initial conditions temperature

% validation temperatures taken from graph 4 final report jasper, from
% excel sheet
t = 0:30:900;
T_thp_measured = [18 28.5 54.25 80 100 120 135 150 162 174 185.5 197 208.5 220 232.5 245 254 263 271 279 278.5 278 276 274 273 272 273.5 275 276.5 278 278];
T_thp_measured_interp = [[0:0.01:900]', interp1(t, T_thp_measured', [0:0.01:900]')];
T_tsensor_measured = [18 28.5 39 60.5 82 102 122 138 154 167 180 192 204 217.5 231 242.5 254 264.5 275 282.5 290 290 290 290 290 287.5 285 288.5 292 290.5 289];

T_bhp_measured = [18 36.95 63 88.5 114 134 154 171 188 203 218 231.5 245 253.5 262 276 290 292 294 290 286 284.5 283 280.5 278 280.5 283 284.5 286 282 278];
T_bhp_measured_interp = [[0:0.01:900]', interp1(t, T_bhp_measured', [0:0.01:900]')];
T_bsensor_measured = [18 28.5 39 59.5 80 102 124 142 160 176 192 206 220 235 250 263 276 282 288 289 290 290 290 290 290 290 290 290 290 290 290];

%% Simulate heat up process
% Simulate heat up process
out = sim('model.slx', 'StartTime','0','StopTime','900','FixedStep','0.01');

%% Plot results
close all
figure;
plot(out.T_thp);
hold on
plot(out.T_tc);
plot(out.T_tm);
plot(out.T_ti);
plot(out.T_r);
title("Top Section Simulated Temperatures")
ylabel("Temperature [C]")
yyaxis right
plot(0:60:900, P_tc_org, '-o');
ylabel("Coil Power [W]")
legend("T thp", "T tc", "T tm", "T ti", "T r", "P tc", 'location', 'e');

figure;
plot(out.T_bhp);
hold on
plot(out.T_bc);
plot(out.T_bm);
plot(out.T_bi);
title("Bottom Section Simulated Temperatures")
ylabel("Temperature [C]")
yyaxis right
plot(0:60:900, P_bc_org, '-o');
ylabel("Coil Power [W]")
legend("T bhp", "T bc", "T bm", "T bi", "P bc", 'location', 'e');

figure;
plot(out.T_thp);
hold on
plot(0:30:900, T_thp_measured, '-o');
plot(0:30:900, T_tsensor_measured, '-o');
ylabel("Temperature [C]")
yyaxis right
plot(0:60:900, P_tc_org, '-o');
ylabel("Coil Power [W]")
legend("T thp simulated", "T thp measured", "T tsensor measured", "Coil Power measured", 'location', 's')
title("Top Section Simulated vs Measured")

figure;
plot(out.T_bhp);
hold on
plot(0:30:900, T_bhp_measured, '-o');
plot(0:30:900, T_bsensor_measured, '-o');
ylabel("Temperature [C]")
yyaxis right
plot(0:60:900, P_bc_org, '-o');
ylabel("Coil Power [W]")
legend("T bhp simulated", "T bhp measured", "T bsensor measured", "Coil Power measured", 'location', 's')
title("Bottom Section Simulated vs Measured")

% Calculate error
%SSE_top = sum((out.T_thp.Data - T_thp_measured_interp(:,2)) .* (out.T_thp.Data - T_thp_measured_interp(:,2)))
%SSE_bottom = sum((out.T_bhp.Data - T_bhp_measured_interp(:,2)) .* (out.T_bhp.Data - T_bhp_measured_interp(:,2)))

%% Minimize SSE by varying h_c
h_c_list = 14:0.25:18;
SSE = [h_c_list', zeros(length(h_c_list), 1)];
for i = 1:length(h_c_list)
    h_c_top = h_c_list(i);
    out = sim('model.slx', 'StartTime','0','StopTime','5400','FixedStep','0.01');
    SSE(i, 2) = sum((out.T_thp.Data - T_thp_measured_interp(:,2)) .* (out.T_thp.Data - T_thp_measured_interp(:,2)));
end

figure;
plot(SSE(:,1), SSE(:,2));
xlabel("h_c [W/(m^2C)]")
ylabel("SSE")
title("SSE as a function of the convective heat transfer coefficient of air")

%% Validation setup for steady state process, no production
% use the heat up cycle from above, then keep a steady input of 150W and
% see how it compares to the measurements
t = 0:60:900;
P_tc_org = [0 600 660 600 660 600 600 660 600 660 360 360 240 300 60 180];
P_tc = [[0:0.01:900]', interp1(t, P_tc_org', [0:0.01:900]')]; % top coil heat added [W]
P_bc_org = [0 660 600 660 600 660 660 600 420 300 240 300 180 240 180 60];
P_bc = [[0:0.01:900]', interp1(t, P_bc_org', [0:0.01:900]')]; % bottom coil head added [W]

t_ss = [900.01:0.01:5400]';
P_tc_ss = [t_ss repmat([150], 450000, 1)];
P_bc_ss = P_tc_ss;
P_tc = [P_tc ; P_tc_ss];
P_bc = [P_bc ; P_bc_ss];

% Set constant 150W power, no initial conditions
T_a = 18;
T_init = T_a;
% P_tc = [0 150; 5400 150];
% P_bc = P_tc;

t = [0:30:900]';
T_thp_measured = [t, [18 28.5 54.25 80 100 120 135 150 162 174 185.5 197 208.5 220 232.5 245 254 263 271 279 278.5 278 276 274 273 272 273.5 275 276.5 278 278]'];
T_bhp_measured = [t, [18 36.95 63 88.5 114 134 154 171 188 203 218 231.5 245 253.5 262 276 290 292 294 290 286 284.5 283 280.5 278 280.5 283 284.5 286 282 278]'];

T_thp_measured = [T_thp_measured; [1800 275; 2700 267; 3600  260; 4500 259; 5400 268]];
T_bhp_measured = [T_bhp_measured; [1800 275; 2700 270; 3600 273; 4500 270; 5400 268]];

T_bhp_measured_interp = [[0:0.01:5400]', interp1(T_bhp_measured(:,1), T_bhp_measured(:,2), [0:0.01:5400]')];
T_thp_measured_interp = [[0:0.01:5400]', interp1(T_thp_measured(:,1), T_bhp_measured(:,2), [0:0.01:5400]')];

%% Simulate steady state process, no production
out = sim('model.slx', 'StartTime','0','StopTime','5400','FixedStep','0.01');

%% Plot results
close all
figure;
plot(out.T_thp);
hold on
plot(out.T_tc);
plot(out.T_tm);
plot(out.T_ti);
plot(out.T_r);
title("Top Section Simulated Temperatures")
ylabel("Temperature [C]")
yyaxis right
plot(P_tc(:,1), P_tc(:,2));
ylabel("Coil Power [W]")
legend("T thp", "T tc", "T tm", "T ti", "T r", "P tc", 'location', 'e');
xlim([0 5400])

figure;
plot(out.T_bhp);
hold on
plot(out.T_bc);
plot(out.T_bm);
plot(out.T_bi);
title("Bottom Section Simulated Temperatures")
ylabel("Temperature [C]")
yyaxis right
plot(P_bc(:,1), P_bc(:,2));
ylabel("Coil Power [W]")
legend("T bhp", "T bc", "T bm", "T bi", "P bc", 'location', 'e');
xlim([0 5400])

figure;
plot(out.T_thp);
hold on
plot(T_thp_measured(:, 1), T_thp_measured(:, 2), '-o');
ylabel("Temperature [C]")
yyaxis right
plot(P_tc(:,1), P_tc(:,2));
ylabel("Coil Power [W]")
legend("T thp simulated", "T thp measured", "Coil Power measured", 'location', 's')
title("Top Section Simulated vs Measured")
xlim([0 5400])

figure;
plot(out.T_bhp);
hold on
plot(T_bhp_measured(:, 1), T_bhp_measured(:, 2), '-o');
ylabel("Temperature [C]")
yyaxis right
plot(P_tc(:,1), P_tc(:,2));
ylabel("Coil Power [W]")
legend("T bhp simulated", "T bhp measured", "Coil Power measured", 'location', 's')
title("Bottom Section Simulated vs Measured")
xlim([0 5400])

figure;
plot(out.T_thp);
hold on
plot(out.T_thp_tf);
ylabel("Temperature [C]")
legend("T bhp differential equations", "T bhp transfer function", 'location', 'south')
title("Comparing transfer function to differential eq model")
%xlim([0 5400])
%ylim([0 300])

figure;
plot(out.T_thp);
hold on
plot(out.T_thp_tf);
ylabel("Temperature [C]")
legend("T thp differential equations", "T thp transfer function", 'location', 'south')
title("Comparing transfer function to differential eq model")
%xlim([0 5400])
%ylim([0 300])

% Calculate error
%SSE_top = sum((out.T_thp.Data - T_thp_measured_interp(:,2)) .* (out.T_thp.Data - T_thp_measured_interp(:,2)))
%SSE_bottom = sum((out.T_bhp.Data - T_bhp_measured_interp(:,2)) .* (out.T_bhp.Data - T_bhp_measured_interp(:,2)))
%SSE_tf = sum((out.T_bhp.Data - out.T_bhp_tf.Data) .* (out.T_bhp.Data - out.T_bhp_tf.Data))

%% estimating transfer function
tt_out = timeseries2timetable(out.T_thp);
tt_out.Properties.VariableNames{1} = 'Output';
tt_in = timeseries2timetable(out.P_tc);
tt_in.Properties.VariableNames{1} = 'Input';
tt = [tt_in tt_out];
sys = tfest(tt, 2);