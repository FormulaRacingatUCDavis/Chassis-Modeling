% Updated 9/25/2020 by Yash Taneja
% Added all SES data for various regions incorporating latest monocoque
% height measurement values, made all comments/formatting consistent, added
% optimization features to cycle through all possible combinations of layup
% schedules (currently quasi-isotropic, but can be easily modified).
%
% Updated 5/22/2020 by Quinn Davies
% Changed some formatting, added section to calculate material properties
% (E, nu, G) for each ply and total laminate.
% Added if statements to control program flow because I didn't need the
% whole script every time.
% Added maximum stress failure criterion.
% Added SES calculations for factor of safety.
%
%
% To-do: Checking work and calculations adding core properties and effect
% of shear stiffness, better failure to work with multiple load cases too.
% Improve work flow and user experience.
%
%
% Classical Lamination Plate Theory Calculator
% Given a laminate's individual plies' material properties and
% thicknesses in their local coordinates, and their orientation relative to
% the global coordinate system, as well as any external mechanical and
% thermal loading, this will calculate the laminate stiffness and local
% stress and strain values.

% Assumes plane strain/stress, so all shear and torsion is in plane, and no
% out of plane properties are considered.

% Refer to Herakovich "Mechanics of Fibrous Composites" Chapters 4 and 5.

% Also, ply and laminae are interchangeable as far as I know.

%% Clean-Up and Setup

close all; clear; clc;

% Figure Interpreter
set(groot, 'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Program Flow Switches
% Because sometimes you don't need all of this and I'm too lazy to make a
% functional UI.

% Plots global stress distributions, needs Switch.Stress.
Switch.Plot = 0;

% Displays outputs.
Switch.Display = 0;

% Calculates local and global stress and strain distributions.
Switch.Stress = 1;

% Calculates load for failure, failure mode and ply, needs Switch.Stress.
Switch.Failure = 1;

% Calculates stuff for SES, needs Switch.Properties.
Switch.SES.SES = 1;

% Calculates engineering constants like elastic modulus for plies/laminate.
Switch.Properties = 1;

%% User Input
% Use this section to input data about materials, laminates, and loading
% Make sure everything is in order and in correct units.

% MATERIAL DATA

% Postscript 1 indicates properties along the fiber or warp direction.
% Postscript 2 indicates properties along the transverse or weft direction.
% Postscript 12 typically indicates an in-plane shear direction, such as
% loading on the 1 side and in the 2 direction, or vice versa.
% Add in material properties in order,  keep the index in mind for
% referencing in the laminate schedule.

Material.Name = {'Gr/Ep','GFRP from MAE237 Project 2', ...
    'BT250E-1_3K8HS368gsm','SG Endgrain Balsa', ...
    'AAA-3.1-1/8-07N-5052 Toray Al HC', 'Toray T700 Uni (P707AG-15)'};

% Average thickness of cured ply [m].
Material.Thickness = [0.25e-3; 0.13e-3; ...
    0.4487e-3; 19.05e-3; ...
    19.05e-3; 0.152e-3];

% Areal density of material [gsm].
Material.Density = [368; 152; ...
    300; 152e3*Material.Thickness(4); ...
    49657.2365*Material.Thickness(5); 1.517e6*Material.Thickness(6)];

% Young's modulus in the fiber direction [Pa].
Material.E1 = [56.25e9; 45e9; 2.31e10; 2.349e9; 5.171e8; 1.28e11];

% Young's Modulus in the transverse direction [Pa].
% Identical to E1 for a balanced fabric (plain weave)
Material.E2 =[56.25e9; 12e9; 2.3856e10; 2.349e9; 5.171e8; 9.03e9];

% Major Poisson's ratio, or the ratio of epsilon2 to epsilon1.
Material.nu12 = [0.31; 0.19; 0.148; 0.3; 0.3; 0.31];

% In-plane shear modulus [Pa].
Material.G12 = [5.8e9; 5.5e9; 3.654e9; 186e6; 4.0334e8; 4.23e9];

% Tensile, compressive, and interplanar shear strengths.
Material.sigma1t = [1; 1020e6; 4.564e8; 14.39e6; 1.448e6*sqrt(3); 2.172e9]; 
Material.sigma2t = [1; 40e6; 4.233e8; 14.39e6; 8.96318e5*sqrt(3); 4.43e7]; 
Material.sigma1c = [1; -620e6; -4.585e8; 10.67e6; 2.068e6; -1.448e9];
Material.sigma2c = [1; -140e6; -4.509e8; 10.67e6; 2.068e6; -2.88e7];
Material.tau12 =   [1; 60e6; 1.0756e8; 2.72e6; 1.172e6; 1.54e8]; % shear

% Thermal expansion coefficient per degrees celsius change.
Material.alpha1 = [0; 3.7e-6; 0; 0; 0; 0];
Material.alpha2 = [0; 30e-6; 0; 0; 0; 0];

% Hygroscopic expansion coefficient per percent moisture content change.
Material.beta1 = [0; 0; 0; 0; 0; 0];
Material.beta2 = [0; 0.2; 0; 0; 0; 0];


% MECHANICAL LOAD DATA

% Input the normal and shear loads per unit length, and bending and torsion
% loads per unit length in the vectors below.  Units should be N/m and Nm/m
% respectively. These can be obtained with FBDs and geometry.

Load.Mech = [0;...% Normal loading in the x direction, N/m.
    0;...% Normal loading in the y direction, N/m.
    0;...% In-plane shear loading, N/m.
    0;...% Bending loading in the x direction, (Nm)/m.
    1;...% Bending loading in the y direction, (Nm)/m.
    0];% In-plane torsion loading, (Nm)/m.


% THERMAL DATA

% Input the expected cure temperature and operating temperature, in
% fahrenheit, or uncomment switch with the commented section for celsius.
Load.Therm.Cure = 29;
Load.Therm.Operating = 20;

% For Fahrenheit Input:
% Load.Therm.Change = 1.8*(Load.Therm.Operating - Load.Therm.Cure - 32);

% For Celsius Input:
Load.Therm.Change = 0; %Load.Therm.Operating - Load.Therm.Cure;


% HYGROSCOPIC DATA

% Input the expected change in % relative humidity from cure to operation.
Load.Hygro.Cure = .29;
Load.Hygro.Operating = .20;
Load.Hygro.Change = 0;
% Load.Hygro.Operating - Load.Hygro.Cure;


% SES DATA

% SES(i).Height: Height of composite panel [m].
% SES(i).EI.Goal: Buckling Modulus Goal (Nm^2).
% SES(i).Strength.Goal: Ultimate Tensile Strength Goal [Pa].
% SES(i).Energy.Goal: Energy absorption goal (assuming brittle/elastic
% 25mm of steel tube bending) [J].
% SES(i).Shear.Goal: Shear Goal (shear calculations inaccurate).
% SES(i).FOS: Factor of Safety Initialization.

% Front Bulkhead Support [FBHS]
SES(1).Height = 274.066/1000;
SES(1).EI.Goal = 4.02e3;
SES(1).Strength.Goal = 9.96e4;
SES(1).Energy.Goal = 1.38e1;
SES(1).Shear.Goal = 0;
SES(1).FOS.EI = 0;
SES(1).FOS.Strength = 0;
SES(1).FOS.Energy = 0;

% Front Hoop Brace [FHB]
SES(2).Height = 273.304/1000;
SES(2).EI.Goal = 1.70e3;
SES(2).Strength.Goal = 4.16e4;
SES(2).Energy.Goal = 5.86;
SES(2).Shear.Goal = 0;
SES(2).FOS.EI = 0;
SES(2).FOS.Strength = 0;
SES(2).FOS.Energy = 0;

% Front Bulkhead [FBH]
SES(3).Height = (min(342.9-264.7, 340.4-264.16))/1000;
SES(3).EI.Goal = 3.40e3;
SES(3).Strength.Goal = 8.32e4;
SES(3).Energy.Goal = 1.17e1;
SES(3).Shear.Goal = 0;
SES(3).FOS.EI = 0;
SES(3).FOS.Strength = 0;
SES(3).FOS.Energy = 0;

% Vertical Side Impact Structure [VSIS]
SES(4).Height = 388.62/1000;
SES(4).EI.Goal = 3.40e3;
SES(4).Strength.Goal = 8.32e4;
SES(4).Energy.Goal = 1.17e1;
SES(4).Shear.Goal = 0;
SES(4).FOS.EI = 0;
SES(4).FOS.Strength = 0;
SES(4).FOS.Energy = 0;

% Horizontal Side Impact Structure [HSIS]
SES(5).Height = 448.056/1000;
SES(5).EI.Goal = 1.70e3;
SES(5).Strength.Goal = 4.16e4;
SES(5).Energy.Goal = 5.86;
SES(5).Shear.Goal = 0;
SES(5).FOS.EI = 0;
SES(5).FOS.Strength = 0;
SES(5).FOS.Energy = 0;

% Main Hoop Brace Support [MHBS]
SES(6).Height = 260.096/1000;
SES(6).EI.Goal = 2.68e3;
SES(6).Strength.Goal = 6.64e4;
SES(6).Energy.Goal = 9.22;
SES(6).Shear.Goal = 0;
SES(6).FOS.EI = 0;
SES(6).FOS.Strength = 0;
SES(6).FOS.Energy = 0;

% Accumulator Side Protection [ASP]
SES(7).Height = 353.8/1000;
SES(7).EI.Goal = 5.11e3;
SES(7).Strength.Goal = 1.25e5;
SES(7).Energy.Goal = 1.76e1;
SES(7).Shear.Goal = 0;
SES(7).FOS.EI = 0;
SES(7).FOS.Strength = 0;
SES(7).FOS.Energy = 0;

% Tractive Side Protection [TSP]
SES(8).Height = 390.398/1000;
SES(8).EI.Goal = 4.02e3;
SES(8).Strength.Goal = 9.96e4;
SES(8).Energy.Goal = 1.38e1;
SES(8).Shear.Goal = 0;
SES(8).FOS.EI = 0;
SES(8).FOS.Strength = 0;
SES(8).FOS.Energy = 0;

% Rear Impact Protection [RIP]
SES(9).Height = 155.448/1000;
SES(9).EI.Goal = 2.68e3;
SES(9).Strength.Goal = 6.64e4;
SES(9).Energy.Goal = 9.22;
SES(9).Shear.Goal = 0;
SES(9).FOS.EI = 0;
SES(9).FOS.Strength = 0;
SES(9).FOS.Energy = 0;


% LAMINATE OPTIMIZATION
Laminate.Optimization.FOS.EI = 4.5;  % Buckling modulus FOS goal.
Laminate.Optimization.FOS.Strength = 4.5;  % Strength FOS goal.
Laminate.Optimization.FOS.Energy = 4.5;  % Energy FOS goal.
Laminate.Optimization.Counter = 1;  % Loop counter initialization.
SES_Optimization = 1:length(SES);  % Tracks unoptimized SES regions.
SES_Optimized = [];  % Tracks optimized SES regions.

% Angles to append to form all combinations (algorithm explained at the end
% of the while loop).
Laminate.Optimization.Seed = [0,45,90];

% Initial ply schedules of one symmetric half, to cycle through. 
% Only quasi-isotropic after these four schedules.
Laminate.Optimization.Array = [0,45,90];

% All ply schedules created (including non-quasi-isotropic).
Laminate.Optimization.ArrayHistory = Laminate.Optimization.Array;

% TempArray is used for intermediate calculations in the algorithm later in
% the script.
Laminate.Optimization.TempArray = [];

% Titles of all outputs collected.
Laminate.Optimization.Outputs = ["", "FBHS", "", "", "FHB", "", "", ...
    "FBH", "", "", "VSIS", "", "", "HSIS", "", "", "MHBS", "", "",...
    "ASP", "", "", "TSP", "", "", "RIP", "", ""; "Ply Schedule",...
    "Buckling Modulus FOS", "Strength FOS", "Energy FOS",...
    "Buckling Modulus FOS", "Strength FOS", "Energy FOS",...
    "Buckling Modulus FOS", "Strength FOS", "Energy FOS",...
    "Buckling Modulus FOS", "Strength FOS", "Energy FOS",...
    "Buckling Modulus FOS", "Strength FOS", "Energy FOS",...
    "Buckling Modulus FOS", "Strength FOS", "Energy FOS",...
    "Buckling Modulus FOS", "Strength FOS", "Energy FOS",...
    "Buckling Modulus FOS", "Strength FOS", "Energy FOS",...
    "Buckling Modulus FOS", "Strength FOS", "Energy FOS"];

start_time = datetime;
disp("Started: " + string(start_time));
disp("Running Optimization...");

while not(isempty(SES_Optimization))  % Until all SES regions optimized.
    % LAMINATE DATA (not a user input)
    % Plies are ordered from bottom to top, or outside to inside.

    % Assigning materials to half of layup schedule.
    % Materials are indexed from 1:length(Material.Name).
    Laminate.Half_Material = [];
    for i = 1:length(Laminate.Optimization.Array(:,1))
        Laminate.Half_Material = [Laminate.Half_Material; 6];
    end

    % References the material for the entire layup schedule.
    Laminate.Material = [ ...
        Laminate.Half_Material; ...
        5; ...  % Core material.
        flip(Laminate.Half_Material) ...
        ];

    % Ply angles relative to the global coordinate system
    Laminate.Angle = [ ...
        Laminate.Optimization.Array(:, Laminate.Optimization.Counter); ...
        0; ...  % Core angle.
        flip(Laminate.Optimization.Array(:,Laminate.Optimization.Counter...
        ))];

    % Number of each specified material/angle.
    % Currently, everything is set at a quantity of 1, but this
    % functionality could be used in the future to allow less work with
    % material/angle specification.
    Laminate.Quantity = [];
    for i = 1:length(Laminate.Angle)
        Laminate.Quantity = [Laminate.Quantity; 1];
    end

    %% Thickness and Density

    N = length(Laminate.Material);  % Number of ply groups in laminate

    % Total thickness of laminate
    Laminate.Thickness.Total = sum(Material.Thickness(Laminate.Material)...
        .* Laminate.Quantity);

    % Z-coordinate of each ply interface, zero being the mid-plane.
    Laminate.z = [0; cumsum(Material.Thickness(Laminate.Material) ...
        .* Laminate.Quantity)] - ones(N+1,1)*Laminate.Thickness.Total/2;

    % Identifies which layer(s) is core based off of thickness of the 
    % material used.
    % Skin layers are very thin, anything thicker than 1 mm is not skin.
    Laminate.Core = Material.Thickness(Laminate.Material) > 1e-3;
    Laminate.Thickness.Core = ...
        sum(Material.Thickness(Laminate.Material(Laminate.Core)) ...
        .* Laminate.Quantity(Laminate.Core));
    Laminate.Thickness.Skin = (Laminate.Thickness.Total ...
        - Laminate.Thickness.Core)/2;

    % Areal density of laminate schedule
    Laminate.Density.Total = sum(Material.Density(Laminate.Material) ...
        .* Laminate.Quantity);
    Laminate.Density.Core = ...
        sum(Material.Density(Laminate.Material(Laminate.Core)) ... 
        .* Laminate.Quantity(Laminate.Core));
    Laminate.Density.Skin = Laminate.Density.Total - Laminate.Density.Core;

    %% Stiffness Matrices
    % Chapter 4
    % This section calculates the elastic properties of the individual 
    % laminae and the total laminate, using the above data.

    % Defines material properties for the each laminae.
    Laminate.Local.E1 = Material.E1(Laminate.Material);
    Laminate.Local.E2 = Material.E2(Laminate.Material);
    Laminate.Local.nu12 = Material.nu12(Laminate.Material);
    Laminate.Local.G12 = Material.G12(Laminate.Material);
    Laminate.Local.alpha1 = Material.alpha1(Laminate.Material);
    Laminate.Local.alpha2 = Material.alpha2(Laminate.Material);
    Laminate.Local.beta1 = Material.beta1(Laminate.Material);
    Laminate.Local.beta2 = Material.beta2(Laminate.Material);

    % Effective mechanical loading due to thermal expansion.
    Load.Therm.Mech = zeros(6,1);

    % Effective mechanical loading due to moisture expansion.
    Load.Hygro.Mech = Load.Therm.Mech;

    % Engineering Strain Conservation
    T.R = [1, 0, 0;...
        0, 1, 0;...
        0, 0, 2];

    % Transformation Matrix 
    % Terms represent the change of each local coordinate direction 
    % relative to the global coordinate directions.
    % Used to rotate material properties, strains, and stresses.
    T.T1 = @(x) [cosd(x)^2, sind(x)^2, 2*sind(x)*cosd(x);...
        sind(x)^2, cosd(x)^2, -2*sind(x)*cosd(x);...
        -sind(x)*cosd(x), sind(x)*cosd(x), cosd(x)^2-sind(x)^2];

    % Uses the funtcion GlobalReducedStiffness. Stiffness is in [Pa].
    [Stiffness.Global, Stiffness.Local, Laminate.Global.alpha, ...
        Laminate.Global.beta, Laminate.Compliance] = ...
        GlobalReducedStiffness(Laminate.Local.E1, Laminate.Local.E2, ...
        Laminate.Local.nu12, Laminate.Local.G12, Laminate.Angle, ...
        Laminate.Local.alpha1, Laminate.Local.alpha2,...
        Laminate.Local.beta1, Laminate.Local.beta2);

    % Laminate normal and shear stiffness matrix A. Normal-Shear coupling
    % terms in A = 0 when off-axis lamina are balanced, i.e. for every
    % theta ply there is a -theta ply, 0 and 90 excepted [Pa*m]/[N/m].
    Stiffness.Normal = zeros(3);

    % Laminate (normal and shear)/(bending and torsion) coupling matrix B
    % B = 0 when laminate is symmetric about midplane [Pa*m^2]/[N].
    Stiffness.Coupling = Stiffness.Normal;

    % Laminate bending and torsion stiffness matrix D
    Stiffness.Bending = Stiffness.Normal;

    % Loops through each ply to sum up laminate stiffness matrices A/B/D,
    % as well as hygrothermal-mechanical loading due to expansion effects
    for i = 1:N
        Stiffness.Normal = Stiffness.Normal + Stiffness.Global(:,:,i) ...
            *(Laminate.z(i+1) - Laminate.z(i));
        Stiffness.Coupling = Stiffness.Coupling+Stiffness.Global(:,:,i) ...
            *(Laminate.z(i+1)^2 - Laminate.z(i)^2)/2;
        Stiffness.Bending = Stiffness.Bending + Stiffness.Global(:,:,i) ...
            *(Laminate.z(i+1)^3 - Laminate.z(i)^3)/3;

        % Equivalent Thermal Force
        Load.Therm.Mech(1:3) = Load.Therm.Mech(1:3) + ...
            Stiffness.Global(:,:,i)*Laminate.Global.alpha(:,i) ...
            *(Laminate.z(i+1) - Laminate.z(i));
        % Equivalent Thermal Moment
        Load.Therm.Mech(4:6) = Load.Therm.Mech(4:6) + ...
            Stiffness.Global(:,:,i)*Laminate.Global.alpha(:,i) ...
            *(Laminate.z(i+1)^2 - Laminate.z(i)^2)/2;
        % Equivalent Hygroscopic Force
        Load.Hygro.Mech(1:3) = Load.Hygro.Mech(1:3) + ...
            Stiffness.Global(:,:,i)*Laminate.Global.beta(:,i) ...
            *(Laminate.z(i+1) - Laminate.z(i));
        % Equivalent Hygroscopic Moment
        Load.Hygro.Mech(4:6) = Load.Hygro.Mech(4:6) + ...
            Stiffness.Global(:,:,i)*Laminate.Global.beta(:,i) ...
            *(Laminate.z(i+1)^2 - Laminate.z(i)^2)/2;
    end

    Stiffness.Laminate = [Stiffness.Normal,Stiffness.Coupling;...
        Stiffness.Coupling,Stiffness.Bending];

    Load.Therm.Mech = Load.Therm.Mech*Load.Therm.Change;
    Load.Hygro.Mech = Load.Hygro.Mech*Load.Hygro.Change;
    Load.Total = Load.Mech + Load.Therm.Mech + Load.Hygro.Mech;

    %% Engineering Properties
    % Calculates the typical engineering properties such as elastic moduli,
    % shear moduli, Poisson's ratios, and other coupling constants for each
    % layer in the laminate in global coordinates.

    if(Switch.Properties)
        C = EngineeringConstants(Laminate.Compliance);
        Laminate.Global.E_x = squeeze(C(1,1,:));
        Laminate.Global.E_y = squeeze(C(2,2,:));
        Laminate.Global.G_xy = squeeze(C(3,3,:));
        Laminate.Global.nu_xy = squeeze(C(1,2,:));
        Laminate.Global.nu_yx = squeeze(C(2,1,:));
        Laminate.Global.eta_x_xy = squeeze(C(1,3,:));
        Laminate.Global.eta_y_xy = squeeze(C(2,3,:));
        Laminate.Global.eta_xy_x = squeeze(C(3,1,:));
        Laminate.Global.eta_xy_y = squeeze(C(3,2,:));

        % Engineering Properties (assuming a symmetric laminate)
        % I think these are pretty pointless, only useful for in plane
        % loading (no bending)
        a = Laminate.Thickness.Total * inv(Stiffness.Normal);
        Laminate.Properties.E_x = 1 / a(1,1);
        Laminate.Properties.E_y = 1 / a(2,2);
        Laminate.Properties.G_xy = 1 / a(3,3);
        Laminate.Properties.nu_xy = -a(2,1) / a(1,1);
        Laminate.Properties.nu_yx = -a(1,2) / a(2,2);
        Laminate.Properties.eta_x_xy = a(1,3) / a(3,3);
        Laminate.Properties.eta_y_xy = a(2,3) / a(3,3);
        Laminate.Properties.eta_xy_x = a(3,1) / a(1,1);
        Laminate.Properties.eta_xy_y = a(3,2) / a(2,2);
    end

    %% Constitutive Relations
    % Chapter 5
    % Takes stiffness and load data to determine in-plane
    % strains and stresses

    if(Switch.Stress)
        % Normal and shear strains at the midplane in (1:3), curvature and
        % twist of the plane in (4:6)
        Strain.Midplane = Stiffness.Laminate\(Load.Total);

        % Strain of each ply in the global coordinate system
        Strain.Global = Strain.Midplane(1:3) + Laminate.z' ...
            .*Strain.Midplane(4:6);

        Strain.Local = zeros(3,2,N);
        Stress.Global = Strain.Local;
        Stress.Local = Strain.Local;

        % Used in the next section to plot stress distributions
        % Stores the bottom and top z-coordinate of each laminate.
        Plot.z = zeros(2*N,1);

        % Finds strains and stresses of each ply in both local and global.
        for i = 1:N
            Strain.Local(:,:,i) = T.R*T.T1(Laminate.Angle(i))/T.R ...
                * Strain.Global(:,i:i+1);
            Stress.Global(:,:,i) = Stiffness.Global(:,:,i) ...
                * (Strain.Global(:,i:i+1) - Load.Therm.Change ...
                * Laminate.Global.alpha(:,i) - Load.Hygro.Change ...
                * Laminate.Global.beta(:,i));
            Stress.Local(:,:,i) = T.T1(Laminate.Angle(i)) ...
                * Stress.Global(:,:,i);
            Plot.z(2*i-1:2*i) = 1000*Laminate.z(i:i+1);
        end

        %% Failure Analysis
        % Uses max stress criterion (aka easiest to compute) to determine 
        % when failure will occur.

        if(Switch.Failure)
            Load.Tension1 = ones(N,1);
            Load.Compression1 = Load.Tension1;
            Load.Tension2 = Load.Tension1;
            Load.Compression2 = Load.Tension1;
            Load.Shear = Load.Tension1;

            % Calculates multiplier of load necessary to induce failure for
            % each ply in the respective failure mode
            for i = 1:N
                Load.Tension1(i) = ...
                    Material.sigma1t(Laminate.Material(i)) ...
                    / max(abs(Stress.Local(1,:,i)));
                Load.Compression1(i) = ...
                    Material.sigma1c(Laminate.Material(i)) ...
                    / max(abs(Stress.Local(1,:,i)));
                Load.Tension2(i) = ...
                    Material.sigma2t(Laminate.Material(i)) ...
                    / max(abs(Stress.Local(2,:,i)));
                Load.Compression2(i) = ...
                    Material.sigma2c(Laminate.Material(i)) ...
                    / max(abs(Stress.Local(2,:,i)));
                Load.Shear(i) = ...
                    Material.tau12(Laminate.Material(i)) ...
                    / max(abs(Stress.Local(3,:,i)));
            end

            Load.Failure = [Load.Tension1,Load.Compression1, ...
                Load.Tension2, Load.Compression2, Load.Shear];
            Load.Critical = Load.Failure(abs(Load.Failure) ...
                == min(abs(Load.Failure),[],'all'));
            Load.Critical = Load.Critical(1);
            [Load.Ply, Load.Mode.Number] = find(Load.Critical ...
                == Load.Failure, 1);

            switch Load.Mode.Number
                case 1
                    Load.Mode.Type = 'Longitudinal Tension';
                case 2
                    Load.Mode.Type = 'Longitudinal Compression';
                case 3
                    Load.Mode.Type = 'Transverse Tension';
                case 4
                    Load.Mode.Type = 'Transverse Compression';
                case 5
                    Load.Mode.Type = 'Shear';
                otherwise
                    Load.Mode.Type = 'Unknown';
            end

            if Switch.Display
                disp(['Failure Mode: ', Load.Mode.Type]);
                disp(['Failure Ply: ', num2str(Load.Ply)]);
                disp(['Failure Load: ', num2str(Load.Critical),' Nm/m']);
            end
        end

        %% Stress Distribution
        % Plots the stress distribution, feel free to add more pretty 
        % pictures and what have you.

        if(Switch.Plot)
            Plot.Stress.Global = zeros(2*N,3);
            Plot.Title.StressGlobal = {'$\sigma_x$','$\sigma_y$', ...
                '$\tau_{xy}$'};
            Plot.Stress.Local = Plot.Stress.Global;
            Plot.Title.StressLocal = {'$\sigma_1$','$\sigma_2$', ...
                '$\tau_{12}$'};

            % Loops through each type of stress (longitudinal, transverse, 
            % shear).
            for j = 1:3
                % Loops through each laminae.
                for i = 1:N
                    Plot.Stress.Global(2*i-1:2*i,j) ...
                        = Stress.Global(j,:,i)/1e6;
                    Plot.Stress.Local(2*i-1:2*i,j) ...
                        = Stress.Local(j,:,i)/1e6;
                end

                figure(1)
                hold on
                subplot(3,1,j)
                plot(Plot.Stress.Global(:,j),Plot.z)
                title(Plot.Title.StressGlobal{j});
                ylabel('z (mm)')
                xlabel('$\sigma$ (MPa)')

                figure(2)
                hold on
                subplot(3,1,j)
                plot(Plot.Stress.Local(:,j),Plot.z)
                title(Plot.Title.StressLocal{j});
                ylabel('z (mm)')
                xlabel('$\sigma$ (MPa)')
            end
        end
    end

    %% SES Calculations
    % Calculates flat plate bending stiffness values for SES reasons

    if(Switch.SES.SES)
        for i = SES_Optimized  % Stops data collection if optimized.
            SES(i).FOS.EI = "";
            SES(i).FOS.Strength = "";
            SES(i).FOS.Energy = "";
        end
        for i = SES_Optimization  % Calculates FOS for unoptimized regions.
            % Buckling Modulus for a laminate, calculated using laminated
            % beam theory from EAE 135 and the Bauchau & Craig book.
            SES(i).EI.Theory = SES(i).Height*sum(Laminate.Global.E_x.*...
                (Laminate.z(2:end).^3 - Laminate.z(1:N).^3))/3;

            % 2nd area moment of inertia of skins, from I = bt^3/12.
            SES(i).I = ((Laminate.Thickness.Core + ...
                Laminate.Thickness.Skin)^3 - Laminate.Thickness.Core^3) ...
                * SES(i).Height/12;

            % Estimated elastic modulus from 3-point bend test.
            SES(i).E = SES(i).EI.Theory / SES(i).I;

            % Factor of safety to reach SES goal of EI.
            SES(i).FOS.EI = SES(i).EI.Theory/SES(i).EI.Goal;

            % Expected slope of load displacement curve of 3-point bend.
            SES(i).Gradient = 48*SES(i).EI.Theory/0.4^3;

            if(Switch.Stress)
                % Load in 3-point bend test required to reach failure
                % moment.
                P = 4*Load.Critical*SES(i).Height/0.4;

                % Ultimate tensile strength of skin calculated using SES
                % method.
                SES(i).Strength.UTS = 0.4*P*Laminate.Thickness.Total ...
                    / SES(i).I/8;

                % Estimated tensile strength of section.
                SES(i).Strength.Theory = 2*Laminate.Thickness.Skin ...
                    * SES(i).Height*SES(i).Strength.UTS;

                % Factor of safety to reach SES goal of strength.
                SES(i).FOS.Strength = SES(i).Strength.Theory ...
                    / SES(i).Strength.Goal;
                SES(i).Energy.Theory = P^2/SES(i).Gradient;
                SES(i).FOS.Energy = SES(i).Energy.Theory ...
                    / SES(i).Energy.Goal;
                SES(i).Shear.Theory = 0.025*pi ...
                    * sum(Material.Thickness(Laminate.Material) ...
                    .* Material.tau12(Laminate.Material));
                SES(i).FOS.Shear = SES(i).Shear.Theory/SES(i).Shear.Goal;

                if Switch.Display
                    disp(['SES strength factor of safety: ', ...
                        num2str(SES(i).FOS.Strength)]);
                    disp(['SES energy factor of safety: ', ...
                        num2str(SES(i).FOS.Energy)]);
                    disp(['SES shear factor of safety: ', ...
                        num2str(SES(i).FOS.Shear)]);
                    disp(['SES buckling modulus factor of safety: ', ...
                        num2str(SES(i).FOS.EI)]);
                end
            end

            % Moves optimized SES regions into SES_Optimized from
            % SES_Optimization.
            % Criteria is FOS > 1 for energy, strength, and buckling mod.
            if SES(i).FOS.Strength ...
                    > Laminate.Optimization.FOS.Strength && ...
                    SES(i).FOS.Energy ...
                    > Laminate.Optimization.FOS.Energy && ...
                    SES(i).FOS.EI ...
                    > Laminate.Optimization.FOS.EI
                SES_Optimization(SES_Optimization == i) = [];
                SES_Optimized = [SES_Optimized, i];
            end
        end
    end

    if Switch.Display
        disp(['Individual skin thickness: ', ...
            num2str(Laminate.Thickness.Skin*1000), ' mm']);
        disp(['Core thickness: ', ...
            num2str(Laminate.Thickness.Core*1000), ' mm']);
        disp(['Laminate thickness: ', ...
            num2str(Laminate.Thickness.Total*1000), ' mm']);
        disp(['Laminate areal density: ', ...
            num2str(Laminate.Density.Total/1000), ' kgsm']);
    end

    % Algorithm to cycle through all possible layup schedules.
    % Counter is intended to count if current layup schedule combinations
    % have been iterated through. If so, new combinations are computed with
    % an added ply.
    if Laminate.Optimization.Counter == size(Laminate.Optimization.Array,2)
        Laminate.Optimization.Counter = 0;  % Counter reset.
        Laminate.Optimization.Array = [];  % Quasi-iso combos array reset.
        
        % Until quasi-isotropic layup schedules are added into
        % Laminate.Optimization.Array (i.e. odd numbers of plies in each
        % half of the layup schedule is not possible).
        while isempty(Laminate.Optimization.Array)
            % Creates four identical copies of each column in
            % Laminate.Optimization.ArrayHistory, and then appends
            % Laminate.Optimization.Seed.
            % This algorithm cycles through all combinations possible.
            for i = 1:length(Laminate.Optimization.ArrayHistory)
                Laminate.Optimization.TempArray = ...
                    cat(2, Laminate.Optimization.TempArray, [ ...
                    Laminate.Optimization.ArrayHistory(:,i), ...
                    Laminate.Optimization.ArrayHistory(:,i), ...
                    Laminate.Optimization.ArrayHistory(:,i); ...
                    Laminate.Optimization.Seed]);
            end

            % Updates Laminate.Optimization.ArrayHistory (includes
            % non-quasi-isotropic layup schedules)
            Laminate.Optimization.ArrayHistory = ...
                Laminate.Optimization.TempArray;

            % Starting with non-quasi-isotropic layup schedules, identify
            % indices to remove and remove them.
            Laminate.Optimization.Array = ...
                Laminate.Optimization.ArrayHistory;
            Laminate.Optimization.TempArray = [];  % Reset TempArray.
            remove_indices = [];  % Indices to remove.

            % In each column, number of 0s and 90s, and number of 45s and
            % -45s must be equal for quasi-isotropic properties.
            for i = 1:length(Laminate.Optimization.Array)
                zeroes = ...
                    length(find(Laminate.Optimization.Array(:,i) == 0));
                nineties = ...
                    length(find(Laminate.Optimization.Array(:,i) == 90));
                fortyfives = ...
                    length(find(Laminate.Optimization.Array(:,i) == 45));

                % Identify indices to remove.
                if ((zeroes + nineties) ~= fortyfives) || (zeroes ~= nineties)
                    remove_indices = [remove_indices, i];
                end
            end

            for i = flip(remove_indices)  % Remove all identified indices.
                Laminate.Optimization.Array(:,i) = [];
            end
            
            % Identify all half-layup schedules that aren't reflected upon
            % their mid-axis. This is an added constraint to avoid negative
            % thermal effects by a two-stage cure at the monocoque's scale.
            remove_indices = [];
            for i = 1:length(Laminate.Optimization.Array)
                trigger = false;
                for j = 1:(size(Laminate.Optimization.Array, 1)/2)
                    if Laminate.Optimization.Array(j,i) ~= ...
                            Laminate.Optimization.Array(end-j+1, i)
                       trigger = true; 
                    end
                end
                if trigger == true
                    remove_indices = [remove_indices, i];
                end
            end
            
            for i = flip(remove_indices)  % Remove all identified indices.
                Laminate.Optimization.Array(:,i) = [];
            end
        end
    end
    clear remove_indices zeroes nineties fortyfives minusfortyfives trigger

    % Store data in Laminate.Optimization.Outputs.
    Laminate.Optimization.Outputs = ...
        [Laminate.Optimization.Outputs; "[" + ...
        strjoin(string(Laminate.Angle(1:length(Laminate.Half_Material)))...
        ,{'/'}) + "/core]s", ...
        SES(1).FOS.EI, SES(1).FOS.Strength, SES(1).FOS.Energy, ...
        SES(2).FOS.EI, SES(2).FOS.Strength, SES(2).FOS.Energy, ...
        SES(3).FOS.EI, SES(3).FOS.Strength, SES(3).FOS.Energy, ...
        SES(4).FOS.EI, SES(4).FOS.Strength, SES(4).FOS.Energy, ...
        SES(5).FOS.EI, SES(5).FOS.Strength, SES(5).FOS.Energy, ...
        SES(6).FOS.EI, SES(6).FOS.Strength, SES(6).FOS.Energy, ...
        SES(7).FOS.EI, SES(7).FOS.Strength, SES(7).FOS.Energy, ...
        SES(8).FOS.EI, SES(8).FOS.Strength, SES(8).FOS.Energy, ...
        SES(9).FOS.EI, SES(9).FOS.Strength, SES(9).FOS.Energy];
    Laminate.Optimization.Counter = Laminate.Optimization.Counter + 1;
end

disp("Completed: " + string(datetime));
disp("Time Elapsed: " + string(datetime - start_time));

%% The Return of the Son of Clean-up's Revenge 2: The Sequelling
% Wouldn't want those index variables cluttering up the workspace,
% now would we?

% Uses Latin Hypercube sampling and conjugate gradient successive over
% relaxation fifth-degree shape function minimization with calculus of
% variations iterated using log(N) operations to solve an elliptic form of
% the Riemann Hypothesis.
clear i j T a C start_time SES_Optimization SES_Optimized

%% Functions

function Q = LocalReducedStiffness(E1,E2,nu12,G12)
% This function calculates the local reduced stiffness matrix Q for an
% orthotropic lamina in plane stress. Input arguments are the four
% independent material properties for this type of material. These
% arguments can be scalars, or vectors of equal length to calculate 
% matrices for an entire laminate.
%
% Significance of input arguments:
% E1: Young's modulus in the fiber direction.
% E2: Young's modulus in the transverse in-plane direction.
% (For a fabric with no bias, these two values are equal).
% nu12: Poisson's ratio relating strain in 2 from stress in 1.
% G12: Shear modulus in plane.
%
% Significance of output arguments:
% Q: Local reduced stiffness matrix, sigma(1,2) = Q*epsilon(1,2).

if ~(isequal(size(E1),size(E2),size(nu12),size(G12)) && isvector(E1) && ...
        isnumeric([E1,E2,nu12,G12]))
    error("Input arguments are not of same size. " + ...
        "Properties must be scalars or numeric vectors.");
end

n = length(E1);

% Allows for storing a 3x3xn array with (:,:,i) being one matrix.
if n > 1
    E1 = reshape(E1,1,1,n);
    E2 = reshape(E2,1,1,n);
    nu12 = reshape(nu12,1,1,n);
    G12 = reshape(G12,1,1,n);
end
nu21 = nu12.*E2./E1;

z = zeros(1,1,n);
denom = 1-nu12.*nu21;

% Local Reduced Stiffness Matrix
Q = [E1, nu21.*E1, z;...
    nu12.*E2, E2, z;...
    z, z, denom.*G12] ./ denom;

end

function [QBar, Q, aBar, bBar, SBar] = ...
    GlobalReducedStiffness(E1, E2, nu12, G12, theta, a1, a2, b1, b2)
% This function calculates the global reduced stiffness matrix Qbar for an
% orthotropic lamina in plane stress.  Also outputs the local reduced
% stiffness matrix Q. Input arguments are the four independent material
% properties for this type of material, as well as the ply angle.  These
% arguments can be scalars, or vectors of equal length to calculate 
% matrices for an entire laminate. If theta is given as [] (no angle
% specified), then it outputs 181 3x3 matrices each for a range of possible
% angles.
%
% Significance of input arguments:
% E1: Young's modulus in the fiber direction.
% E2: Young's modulus in the transverse in-plane direction (For a fabric
% with no bias, these two values are equal).
% nu12: Poisson's ratio relating strain in 2 from stress in 1.
% G12: Shear modulus in plane.
% theta: Angle between local, fiber direction 1 and global, structure
% direction x.
%
% Significance of output arguments:
% Q: Local reduced stiffness matrix, sigma(1,2) = Q*epsilon(1,2).
% QBar: Global reduced stiffness matrix, sigma(x,y) = QBar*epsilon(x,y).
% SBar: Global reduced compliance matrix, epsilon(x,y) = SBar*sigma(x,y).

if ~(isequal(size(E1),size(E2),size(nu12),size(G12)) && isvector(E1) && ...
        isnumeric([E1,E2,nu12,G12]) && (isequal(size(E1),size(theta)) ||...
        (isempty(theta) && length(E1) == 1)))
    error("Input arguments are not of same size. " + ...
        "Properties must be scalars or numeric vectors.");
end

% Local Reduced Stiffness Matrix
Q = LocalReducedStiffness(E1,E2,nu12,G12);

% If no angle given, calculates QBar over range of possible orientations.
if isempty(theta)
    theta = reshape(-90:90,1,1,181);
    Q = Q.*ones(1,1,181);
    N = 181;
else
    N = length(theta);
    theta = reshape(theta,1,1,N);
end

% Engineering Strain Conservation
R = [1, 0, 0;...
    0, 1, 0;...
    0, 0, 2];

% Tation Matrix
% Terms represent the change of each local coordinate direction relative
% to the global coordinate directions.
T1 = @(x) [cosd(x)^2, sind(x)^2, 2*sind(x)*cosd(x);...
    sind(x)^2, cosd(x)^2, -2*sind(x)*cosd(x);...
    -sind(x)*cosd(x), sind(x)*cosd(x), cosd(x)^2-sind(x)^2];

% Initializes array.
QBar = zeros(3,3,N);
aBar = zeros(3,N);
bBar = aBar;
SBar = QBar;

for i = 1:N
    % Global Reduced Stiffness Matrix
    QBar(:,:,i) = T1(-theta(i))*Q(:,:,i)*R*T1(theta(i))/R;
    aBar(:,i) = R*T1(-theta(i))*[a1(i);a2(i);0];
    bBar(:,i) = R*T1(-theta(i))*[b1(i);b2(i);0];
    SBar(:,:,i) = inv(QBar(:,:,i));
end
end

function [C] = EngineeringConstants(S)
% This function calculates the engineering material properties for an
% orthotropic lamina in plane stress.

% Significance of input arguments:
% S: Global reduced compliance matrix, epsilon(x,y) = S*sigma(x,y).

% C: Array of global material properties, which behave like the input ones,
% but in the global coordinate system instead of the local coordinates.
% Not present in local coordinates are the normal shear coupling terms,
% eta_{i,ij} and eta_{ij,i}.  These relate shear strains to normal stresses
% and normal strains to shear stresses.
%
% Equations for all functions taken from Carl T. Herakovich
% "Mechanics of Fibrous Composites" Chapter 4.

% Global Lamina Engineering Constants
Ex = S(1,1,:).^-1;
Ey = S(2,2,:).^-1;
Gxy = S(3,3,:).^-1;
nuxy = -S(1,2,:)./S(1,1,:);
nuyx = -S(2,1,:)./S(2,2,:);

% Coefficients of Mutual Influence (Normal-Shear Coupling Terms)
% First index group is induced strain.
% Second index group is applied stress.
% For example, x_xy would relate normal strain induced by an applied shear
% stress.
etax_xy = S(1,3,:)./S(3,3,:);
etaxy_x = S(3,1,:)./S(1,1,:);
etay_xy = S(2,3,:)./S(3,3,:);
etaxy_y = S(3,2,:)./S(2,2,:);

C = [Ex, nuxy, etax_xy;...
    nuyx, Ey, etay_xy;
    etaxy_x, etaxy_y, Gxy];

end
