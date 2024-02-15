  clc; clear; close all;  

%% Vehicle Parameters 

ms = 50;
mu = 10;
mt = 4;
ks = 50000; 
kt = 127000; 
d = 0.5;   

%% Run Model / Designat Model Input 

Out = sim('QuarterCar.slx');