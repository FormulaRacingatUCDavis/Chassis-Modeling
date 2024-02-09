clc
clear
close all
%% Quarter Car Modelling - Speed Hump Profile
% This generates a speed hump profile as an input signal into the quarter
% car model to estimate suspension and chassis response to the speed hump
% profile.
% 
% Inputs:
%   H - height of the hump in cross-section view [m]
%   L - is the width of the hump in the cross-section view [m]
%   v - is the velocity of the car passing speed hump [m/s]
%   t - is the elapsed time that the car passed the speed hump [s]
% 
% Outputs:
%
% Author(s);
% Tristan Pham       (atlpham@ucdavis.edu)
%
% Last Updated: 6-Sept-2022

%%% Parameters
H = 0.075; %0.3556;
L = 0.5;   %0.2667;
v = linspace(0,2,1000);
t = linspace(0,1,1000);
x = v.*t;
z = linspace(0,1,1000);
z(end) = 0;

%%% Counters
m = 1;
n = 1;
b = 1;
c = 1;

%%% Plotting Equation
while t(m) ~= t(end)
    if (x(n) >= 0) && (x(n) <= L)
        z(b) = (-0.5)*H;
        z(b) = z(b)*(cos(2*pi*((v(c)*t(m))/L))-1);
        m = m + 1;
        n = n + 1;
        b = b + 1;
        c = c + 1;
    else 
        z(b) = 0;
        m = m + 1;
        n = n + 1;
        b = b + 1;
        c = c + 1;
    end
end

%%% Plotting
plot(t,z)
