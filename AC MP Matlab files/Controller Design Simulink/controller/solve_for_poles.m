clear all 
close all
clc

%solve for poles in system equation

coef_vector=[1,0.0288,(2.3163*10^-4),(5.8211*10^-7),(2.4497*10^-10)];
x=roots(coef_vector)

H= tf([0.0030,(4.8978*10^-5),(2.5849*10^-7),(4.3592*10^-10)],[1,0.0288,(2.3163*10^-4),(5.8211*10^-7),(2.4497*10^-10)])
rlocus(H)
sgrid

%calculate damping ratio from max OS 15%
OS=45;
t_settling=250;

damping_ratio= (-log(OS/100))/(sqrt((pi^2)+(log((OS/100)^2))))
omega_n= - ((log(0.02))/(damping_ratio*t_settling))
desired_pole1= -(omega_n*damping_ratio)+(omega_n*sqrt(1-(damping_ratio^2)))*i
desired_pole2= -(omega_n*damping_ratio)-(omega_n*sqrt(1-(damping_ratio^2)))*i

%pole real axis location
P_re= omega_n*damping_ratio;
P_i=(omega_n*sqrt(1-(damping_ratio^2)))
%root locus angle condition

theta_1= 180- (atan((P_re-0.0005)/P_i)*(180/pi))
theta_2= 180- (atan((P_re-0.0038)/P_i)*(180/pi))
theta_3= 180- (atan((P_re-0.0072)/P_i)*(180/pi))
theta_4= (atan((0.0173-P_re)/P_i)*(180/pi))

