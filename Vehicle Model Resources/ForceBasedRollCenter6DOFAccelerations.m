<<<<<<< Updated upstream
function [LongAcc, LatAcc, YawAcc, LongAccTot, LatAccTot, VertAcc, RollAcc, PitchAcc] = ...
    ForceBasedRollCenter6DOFAccelerations( TFx, TFy, TMz, AFx, AFy, AMy, AMz, AFz, ...       % Loads
        Wheelbase, TrackWidth, TrackWidthTest, Steer, PercentFront, ...                      % Geometry
        Front_View_Suspension_Arm, z, Side_View_Suspension_Arm, Front_View_Jacking_Angle, ...
        Side_View_Jacking_Angle, ...
        Front_Roll_Stiffness, Rear_Roll_Stiffness, Roll_Dampening_Coeff_Front, ...
        Roll_Dampening_Coeff_Rear, Roll_Speed, Roll, ...
        Pitch_Stiffness, Pitch_Dampening, Pitch, Pitch_Velocity, ...
        Ride_Stiffness, Ride, Ride_Dampening, Ride_Velocity,...
        Mass, YawInertia, CoG, Sprung_COG_Height, ...                                       % Inertia
        LongVel, LatVel, YawVel, Ixx, Iyy)                                                  % Velocities
%% FroceBasedRollCenter6DOFAccelerations - Force Based Roll Center Accelerations
% Computes force based roll center accelerations for a 6DOF full track chassis model
% https://www.overleaf.com/project/5e73ad44ea180c00018cb332 (Section ?)
%
=======
function [LongAcc, LatAcc, YawAcc, LongAccTot, LatAccTot, VertAcc] = ...
    ForceBasedRollCenter6DOFAccelerations( TFx, TFy, TMz, AFx, AFy, AMz, ... % Loads
        Wheelbase, TrackWidth, Steer, ...                         % Geometry
        Mass, YawInertia, CoG, ...                                % Inertia
        LongVel, LatVel, YawVel )                                 % Velocities
%% FroceBasedRollCenter6DOFAccelerations - Force Based Roll Center Accelerations
% Computes planar motion accelerations for a 3DOF full track chassis model
% 
>>>>>>> Stashed changes
% Inputs:
%   TFx        - (n,4 numeric) Tire Longitudinal Force {T_F_y}    [N]
%   TFy        - (n,4 numeric) Tire Lateral Force      {T_F_y}    [N]
%   TMz        - (n,4 numeric) Tire Aligning Moment    {T_M_z}    [N-m]
%   AFx        - (n,1 numeric) Aerodynamic Drag Force  {A_F_x}    [N]
%   AFy        - (n,1 numeric) Aerodynamic Side Force  {A_F_y}    [N]
%   AMz        - (n,1 numeric) Aerodynamic Yaw Moment  {A_M_z}    [N-m]
%   Wheelbase  - (n,1 numeric) Wheelbase               {L}        [m]
%   TrackWidth - (n,2 numeric) Track Width             {t_w}      [m]
%   Steer      - (n,4 numeric) Tire Steer Angle        {delta}    [deg]
%   Mass       - (n,1 numeric) Total Vehicle Mass      {m}        [kg]
%   YawInertia - (n,1 numeric) Vehicle Yaw Inertia     {I_zz}     [kg-m^2]
%   CoG        - (n,3 numeric) Center of Gravity       {CoG}      [m]
%   LongVel    - (n,1 numeric) Longitudinal Velocity   {dot{x}}   [m/s]
%   LatVel     - (n,1 numeric) Lateral Velocity        {dot{y}}   [m/s]
%   YawVel     - (n,1 numeric) Yaw Velocity            {dot{psi}} [rad/s]
%
% Outputs:
%   LongAcc    - (n,1 numeric) Longitudinal Acceleration       {ddot{x}}  [m/s^2]
%   LatAcc     - (n,1 numeric) Lateral Acceleration            {ddot{y}}  [m/s^2]
%   YawAcc     - (n,1 numeric) Yaw Acceleration                {ddot{psi} [rad/s^2]
%   LongAccTot - (n,1 numeric) Total Longitudinal Acceleration {a_x}      [m/s^2]
%   LatAccTot  - (n,1 numeric) Total Lateral Acceleration      {a_y}      [m/s^2]
%
% Notes:
<<<<<<< Updated upstream
% check for correct yaw acc implementation
% Itterations need to be checked for rotational accelerations, like tire
% forces and yaw moments at each tire
% How do we solve for all suspension parameters, and how do they change
% over time
% Test Case fix
% Vertical Acc Parameters fix, understanding ride and implementing it
% IXX and IYY from CAD
% Array operations
=======
>>>>>>> Stashed changes
%
% Author(s): 
% Tristan Pham       (atlpham@ucdavis.edu)        [Oct 2020 - ???     ]
% 
<<<<<<< Updated upstream
% Last Updated: 7-May-2022
=======
% Last Updated: 17-Feburary-2022
>>>>>>> Stashed changes

%% Test Case
if nargin == 0
    TFx = zeros(1,4);
    TFy = [250, 250, 0, 0];
    TMz = zeros(1,4);
    
    AFx = -100;
    AFy = 0;
    AMz = 0;
<<<<<<< Updated upstream
    AFz = 0;
    AMy = 0;
    
    Wheelbase  = 1.575; 
    TrackWidth = 1.22;               %*ones(1,2); (FIXME)did a quick fix for array sizes 
    TrackWidthTest = 1.22*ones(1,2);
    Steer      = [18, 15, 1, -1];
    
    Mass         = 275;   
    YawInertia   = 130; 
    CoG          = [(0.47-0.5)*Wheelbase, 0, 0.25];
    PercentFront = 0.5;
    Sprung_COG_Height = 0;
=======
    
    Wheelbase  = 1.575; 
    TrackWidth = 1.22*ones(1,2);
    Steer      = [18, 15, 1, -1];
    
    Mass       = 275;   
    YawInertia = 130; 
    CoG        = [(0.47-0.5)*Wheelbase, 0, 0.25];
>>>>>>> Stashed changes
    
    LongVel = 15;
    LatVel  = 0;
    YawVel  = 0.5;
<<<<<<< Updated upstream

    Front_View_Suspension_Arm = 0;   %calculation
    Side_View_Suspension_Arm = 0;
    Front_View_Jacking_Angle = 0;    
    Side_View_Jacking_Angle = 0;    
    
    Front_Roll_Stiffness = 0;
    Rear_Roll_Stiffness = 0;
    Roll_Dampening_Coeff_Front = 0;
    Roll_Dampening_Coeff_Rear = 0;
    Roll_Speed = 0;
    Roll = 0;

    Pitch_Stiffness = 0;
    Pitch_Dampening = 0;
    Pitch = 0;
    Pitch_Velocity = 0;

    Ride_Stiffness = 0;
    Ride = 0;
    Ride_Dampening = 0;
    Ride_Velocity = 0;

    Iyy = 0;
    Ixx = 0;

    z = 0;
    
    [LongAcc, LatAcc, YawAcc, LongAccTot, LatAccTot, VertAcc, RollAcc, PitchAcc] = ...
    ForceBasedRollCenter6DOFAccelerations( TFx, TFy, TMz, AFx, AFy, AMy, AMz, AFz, ...      
        Wheelbase, TrackWidth, TrackWidthTest, Steer, PercentFront, ...                      
        Front_View_Suspension_Arm, z, Side_View_Suspension_Arm, Front_View_Jacking_Angle, ...
        Side_View_Jacking_Angle, ...
        Front_Roll_Stiffness, Rear_Roll_Stiffness, Roll_Dampening_Coeff_Front, ...
        Roll_Dampening_Coeff_Rear, Roll_Speed, Roll, ...
        Pitch_Stiffness, Pitch_Dampening, Pitch, Pitch_Velocity, ...
        Ride_Stiffness, Ride, Ride_Dampening, Ride_Velocity,...
        Mass, YawInertia, CoG, Sprung_COG_Height, ...                                      
        LongVel, LatVel, YawVel, Ixx, Iyy)                    
=======
    
    [LongAcc, LatAcc, YawAcc, LongAccTot, LatAccTot] = ...
        ForceBasedRollCenter6DOFAccelerations( ...
            TFx, TFy, TMz, AFx, AFy, AMz, ...
            Wheelbase, TrackWidth, Steer, ... 
            Mass, YawInertia, CoG, ... 
            LongVel, LatVel, YawVel ) %#ok<NOPRT>
>>>>>>> Stashed changes
    
    return;
end

%% Computations

<<<<<<< Updated upstream
%%% Tire Positions (n,4,3 numeric)
TirePos = cat( 3, Wheelbase/2 + [-1, -1, 1, 1].*CoG(:,1), ...
          [TrackWidthTest(:,1).*[1 -1]/2, TrackWidthTest(:,2).*[1 -1]/2], ...
          zeros( size(Steer) ) );
=======
%%% Matricies (FIXME)
x_vector = [x;y;z];
omega_vector = [roll,yaw,pitch];
Is = [Ixx_s, 0, -Ixz_s; 0, Iyy_s, 0; -Izx_s, 0, Izz_s];

[sum_Fx;sum_Fy;sum_Fz] = ms.*(x_vector_ddot + cross(omega_vector_dot,x_vector_dot)) + (muf + mur) .* ...
    ([x_ddot;y_ddot;0] + cross([0;0;psi_dot],[x_dot;y_dot;0]));
    
[sum_Mx;sum_My;sum_Mz] = Is .* omega_vector_ddot + cross(omega_vector_dot, Is .* Omega_vector_dot) ...
    + [0;0;Izz_s.*psi_ddot];

mu(i)*zu_ddot(i) = kt(i).*(zr(i) - zu(i)) - kr(i).*(zu(i)-zc(i)) - ...
    br(i).*(zu_dot(i)-zc_dot(i));

%%% Tire Positions (n,4,3 numeric)
TirePos = cat( 3, Wheelbase/2 + [-1, -1, 1, 1].*CoG(:,1), ...
                  [TrackWidth(:,1).*[1 -1]/2, TrackWidth(:,2).*[1 -1]/2], ...
                  zeros( size(Steer) ) );
>>>>>>> Stashed changes

%%% Tire Loads (n,4,3 numeric)
TireLoad = cat( 3, TFx, TFy, zeros( size(Steer) ) );

%%% Tire Yaw Moment (n,1 numeric)
TireMoment = cross(TirePos, TireLoad, 3);
<<<<<<< Updated upstream
TMz = sum( TireMoment(:,:,3) ) + sum(TMz,2);

%%%Wheelbase Calculations
a = Wheelbase .* (1-PercentFront);
b = Wheelbase .* PercentFront;

%%%Chassis Accelerations

%%Yaw Acceleration
%Note: in overleaf there is the option of a or b

YawComputation1 = sum((TFx .* cosd(Steer) - TFy .* sind(Steer)) .* -TrackWidth ./ 2 ...
    + (TFx .* sind(Steer) + TFy .* cosd(Steer)) .* (a) + TMz);

YawComputation2 = sum((TFx .* cosd(Steer) - TFy .* sind(Steer)) .* TrackWidth ./ 2 ...
    + (TFx .* sind(Steer) + TFy .* cosd(Steer)) .* (a) + TMz);

YawAcc = ((YawComputation1 + YawComputation2) - AFy .* (Wheelbase ./ 2 - a) - AMz) ./ YawInertia;

%%Roll Acceleration 
% What to put for z
RollAcc = (sum((TFx .* sind(Steer) + TFy .* cosd(Steer)) .* Front_View_Suspension_Arm...
    - (TFx .* cosd(Steer) - TFy .* sind(Steer)) .* tand(Side_View_Jacking_Angle) .* TrackWidth ./ 2)...
    - AFy .* (Sprung_COG_Height + z) - AMz - (Front_Roll_Stiffness + Rear_Roll_Stiffness) .* Roll...
    - (Roll_Dampening_Coeff_Front + Roll_Dampening_Coeff_Rear) .* Roll_Speed) ./ Ixx;

%%Pitch Acceleration
%Note* in overleaf there is the option of a or b
PitchAcc = (sum(-(TFx .* cosd(Steer) - TFy .* sind(Steer)) .* Side_View_Suspension_Arm...
    + (TFx .* sind(Steer) + TFy .* cosd(Steer)) .* tand(Front_View_Jacking_Angle).* (a) )...
    + AFx .* (Sprung_COG_Height + z) - AFz .* (a - Wheelbase ./ 2) - AMy ...
    - Pitch_Stiffness .* Pitch - Pitch_Dampening .* Pitch_Velocity) ./ Iyy;


%%Longitudinal Acceleration
LongAcc = ((sum( TFx.*cosd(Steer) - TFy.*sind(Steer), 2 ) - AFx) ./ Mass) + LatVel.*YawVel;

%%Lateral Acceleration
LatAcc  = ((sum( TFy.*cosd(Steer) + TFx.*sind(Steer), 2 ) - AFy) ./ Mass) - LongVel.*YawVel;

%%Vertical Acceleration
VertAcc = (sum(Ride_Stiffness .* Ride + Ride_Dampening .* Ride_Velocity + (TFx .* cosd(Steer) - TFy .* sind(Steer)) .* tand(Side_View_Jacking_Angle)...
    + (TFy .* cosd(Steer) + TFx .* sind(Steer)) .* tan(Front_View_Jacking_Angle)) - AFz) ./ Mass;

LongAccTot = LongAcc - LatVel .*YawVel;
LatAccTot  = LatAcc  + LongVel.*YawVel;














%%% Matricies (FIXME) Other way of expressing computations
% x_vector = [x;y;z];
% omega_vector = [roll,yaw,pitch];
% Is = [Ixx_s, 0, -Ixz_s; 0, Iyy_s, 0; -Izx_s, 0, Izz_s];
% 
% sum_F = ms.*(x_vector_ddot + cross(omega_vector_dot,x_vector_dot)) + (muf + mur) .* ...
%     ([x_ddot;y_ddot;0] + cross([0;0;psi_dot],[x_dot;y_dot;0]));
%     
% sum_M = Is .* omega_vector_ddot + cross(omega_vector_dot, Is .* Omega_vector_dot) ...
%     + [0;0;Izz_s.*psi_ddot];
% 
% mu(i)*zu_ddot(i) = kt(i).*(zr(i) - zu(i)) - kr(i).*(zu(i)-zc(i)) - ...
%     br(i).*(zu_dot(i)-zc_dot(i));
=======
TireMoment = sum( TireMoment(:,:,3) ) + sum(TMz,2);

%%% Moments
%Correct a or b
Mx = sum((TFx.*sind(Steer) + TFy.*cosd(Steer)).*Front_View_Suspension_Arm...
    - (TFx.*cosd(Steer)-TFy.*sind(Steer)).*tand(Side_View_Jacking_Angle).*tw/2)...
    - AFy .* (Sprung_COG_Height + z) - AMz - (front_roll_stiffness + rear_roll_stiffness).*roll...
    - (roll_dampening_coeff_front + roll_dampening_coeff_rear).*roll_speed;

%Correct a or b
My = sum(-(TFx.*cosd(Steer) - TFy.*sind(Steer)).*Side_View_Suspension_Arm...
    + (TFx.*sind(Steer) + TFy.*cosd(Steer)).*tand(Front_View_Jacking_Angle).* (a or b) )...
    + AFx .* (Sprung_COG_Height + z) - AFz.*(a-L/2) - AMy - pitch_stiffness.*pitch - pitch_dampening.*pitch_velocity;

%Positive or negative trackwidth?
Mz = sum((TFx.*cosd(Steer) - TFy.*sind(Steer)).*(-)tw/2 ...
    + (TFx.*sind(Steer) + TFy.*cosd(Steer)).*(a or -b) +TMz)...
    - AFy .* (L/2 - a) - AMz; 

%%% Chassis Accelerations
LongAcc = ((sum( TFx.*cosd(Steer) - TFy.*sind(Steer), 2 ) - AFx) ./ Mass) + LatVel.*YawVel;

LatAcc  = ((sum( TFy.*cosd(Steer) + TFx.*sind(Steer), 2 ) - AFy) ./ Mass) - LongVel.*YawVel;

VertAcc = (sum(kr + br.*j + (TFx.*cosd(Steer) - TFy.*sind(Steer).*tand(Side_View_Jacking_Angle)...
    + (TFy.*cosd(Steer) + TFx.*sind(Steer)).*tan(Front_View_Jacking_Angle) - AFz) ./ Mass;

YawAcc  = (TireMoment - AMz) ./ YawInertia;

LongAccTot = LongAcc - LatVel .*YawVel;
LatAccTot  = LatAcc  + LongVel.*YawVel;
>>>>>>> Stashed changes
