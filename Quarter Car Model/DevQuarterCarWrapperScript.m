  clc; clear; close all;  

%% Vehicle Parameters 

Parameter.SuspensionSpringStiff = 80000;
Parameter.TireStiff             = 87.45;

Parameter.SprungMass   = 2500;
Parameter.UnsprungMass = 320;

Parameter.DampingCoeff = 350;

%% Run Model / Designat Model Input 

Out = sim('QuarterCar.slx');