clc; clear; close all;

% Impact Attenuator Optimization Script
% Blake Christierson
% bechristierson@ucdavis.edu
% 6-17-2020

% Figure Interpreter
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Constraint Specs
Spec.SF = 1.5;

Spec.AvgAcc = 20; % Maximum Average Deceleration [g]
Spec.MaxAcc = 40; % Maximum Peak Deceleration [g]

Spec.Energy = 7350; % Minimum Absorbed Energy [J]

Spec.Mass = 300; % Mass of Vehicle in Impact Case [kg]
Spec.Velocity = 7; % Velocity of Vehicle in Impact Case [m/s]

Spec.Dim = [100 100 100  ; ...
            150 175 200]; % [Width, Height, Thickness] Dimension Bounds [mm]

% Conversion to MKS
Spec.AvgAcc = Spec.AvgAcc .* 9.81; % Maximum Average Deceleration [m/s^2]
Spec.MaxAcc = Spec.MaxAcc .* 9.81; % Maximum Peak Deceleration [m/s^2]

Spec.Dim = Spec.Dim ./ 1000; % [Width, Height, Thickness] Dimension Bounds [m]

%% Material Properties
CrushLite.Rho = [ 0.6 1   1.2 1   1, ...
                  1.6 1.8 1.6 2   2.3, ...
                  2.3 2.3 3   3.6 3.1, ...
                  3.4 3.1 3.7 4.2 4.3, ...
                  5.2 4.5 4.5 5.2 5.4, ...
                  5.7 6.0 5.7 6.1 6.1, ...
                  8.1 8.1]'; % Density [lb/ft^3]
              
CrushLite.CS = [  7.5 10  25  25  35, ...
                  45  45  50  75  80, ...
                  90  100 120 120 120, ...
                  140 170 180 210 230, ...
                  245 275 320 330 350, ...
                  380 420 440 450 535, ...
                  700 750 ]'; % Crush Strength [psi]
              
CrushLite.MS = ones( size( CrushLite.CS ) ) .* 0.7; % Minimum Strain [m/m]

% Conversion to MKS
CrushLite.Rho = CrushLite.Rho .* 16.0185; % Density [kg/m^3]
CrushLite.CS = CrushLite.CS .* 6.89476; % Crush Strengh [kPa]

%% Declare Optimization Structures
Solution = struct();

%% Loop Through Optimization Problem For Varying Materials
for i = 1 : length( CrushLite.Rho )
    
%% Declare Optimization Variables
Variable.Width = optimvar( 'Width', 'LowerBound', Spec.Dim(1,1), ...
    'UpperBound', Spec.Dim(2,1) );
Variable.Height = optimvar( 'Height', 'LowerBound', Spec.Dim(1,2), ...
    'UpperBound', Spec.Dim(2,2) );
Variable.Thickness = optimvar( 'Thickness', 'LowerBound', Spec.Dim(1,3), ...
    'UpperBound', Spec.Dim(2,3) );

%% Declare Constraints & Inequalities
Constraint.Energy = Spec.Energy <= Spec.Mass .* Spec.AvgAcc ./ Spec.SF .* ...
    Variable.Thickness .* CrushLite.MS(i); % Minimum Energy Constraint

Constraint.Force = Spec.Mass .* Spec.AvgAcc ./ Spec.SF >= ...
    CrushLite.CS(i) .* 1000 .* Variable.Width .* Variable.Height; % Crush Strength Constraint

%% Declare Optimization Problem
Problem = optimproblem( 'ObjectiveSense', 'minimize');

Problem.Objective = Variable.Width .* Variable.Height .* ...
    Variable.Thickness .* CrushLite.Rho(i); % Minimizing Weight [kg]

Problem.Constraints.Energy = Constraint.Energy;
Problem.Constraints.Force = Constraint.Force;

%% Finding Initial Condition
x0.Width = sqrt( 2*Spec.Mass*Spec.AvgAcc/(CrushLite.CS(i)*1000*Spec.SF) );
x0.Height = x0.Width / 2;
x0.Thickness = Spec.Energy .* Spec.SF ./ ...
    (Spec.Mass .* Spec.AvgAcc .* CrushLite.MS(i) );
    
%% Solving Optimization Problem
[Sol, Weight] = solve( Problem, x0 );

Sol.Weight = Weight;
if isempty( fieldnames(Solution) )
    Solution = Sol;
else
    Solution(i) = Sol;
end

end

clear i Variable Constraint Problem x0 Weight

%% Plotting
figure
subplot(4,1,1)
plot( CrushLite.CS, [Solution.Height] .* 39.37 )
hold on
plot( CrushLite.CS, [Solution.Width] .* 39.37 )
plot( CrushLite.CS, [Solution.Thickness] .* 39.37 )

xlabel( 'Crushlite Crush Strength [$kPa$]' )
ylabel( 'Dimension [$in$]' )
legend( 'Height', 'Width', 'Thickness' )

subplot(4,1,2)
plot( CrushLite.CS, [Solution.Weight] )

xlabel( 'CrushLite Crush Strength [$kPa$]' )
ylabel( 'Weight [$kg$]' )

subplot(4,1,3)
plot( CrushLite.CS, Spec.Energy .* ones( size(CrushLite.CS) ), 'k--' )
hold on
plot( CrushLite.CS, Spec.Mass .* Spec.AvgAcc .* ...
    [Solution.Thickness]' .* CrushLite.MS )

xlabel( 'CrushLite Crush Strength [$kPa$]' )
ylabel( 'Energy Absorption [$J$]' )
legend( 'Minimum Energy', 'Solution Energies' )

subplot(4,1,4)
plot( CrushLite.CS, Spec.AvgAcc ./ 9.81 .* ones( size(CrushLite.CS) ), 'k--' )
hold on
plot( CrushLite.CS, CrushLite.CS .* 1000 .* ...
    [Solution.Width]' .* [Solution.Height]' ./ (Spec.Mass .* 9.81 ) )

xlabel( 'CrushLite Crush Strength [$kPa$]' )
ylabel( 'Average Deceleration [$g$]' )
legend( 'Max Average Decel', 'Solution Decels' )