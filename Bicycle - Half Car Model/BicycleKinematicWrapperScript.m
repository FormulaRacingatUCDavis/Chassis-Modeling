  clc; clear; close all;  

%% Vehicle Parameters 

L = 1.525;
I   = 2.75;
m   = 275;
PercentFront = 0.6;
k = 400;

%% Run Model / Designat Model Input 

Out = sim('SecondKinematic.slx');