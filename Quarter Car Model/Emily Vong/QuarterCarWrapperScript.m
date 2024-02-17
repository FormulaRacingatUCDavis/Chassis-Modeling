clc; clear; close all;

m_s = 50; % sprung mass in kg (mass of the car body)
m_u = 10; % unsprung mass in kg (mass of the car body)
m_t = 4; % tire mass in kg
k_s = 50000; % spring constant of suspension in N/m
k_t = 127000; % spring constant of tire in N/m
d = 50; %damping coefficient in Ns/m
h_s = 1; %0.762 % height of sprung mass
h_u = 0.5; %0.1524 % height of unsprung mass
h_t = 0; %0.1 % height of tire

Out = sim("QuarterCarModel1.slx")
