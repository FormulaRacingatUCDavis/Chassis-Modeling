function [LongAcc, LatAcc, YawAcc, LongAccTot, LatAccTot, VertAcc, RollAcc, PitchAcc] = ...
    ForceBasedRollCenter6DOFAccelerations( TFx, TFy, TMz, AFx, AFy, AMz, ... % Loads
        Wheelbase, TrackWidth, Steer, ...                         % Geometry
        Mass, YawInertia, CoG, ...                                % Inertia
        LongVel, LatVel, YawVel)                                 % Velocities
%% FroceBasedRollCenter6DOFAccelerations - Force Based Roll Center Accelerations
% Computes planar motion accelerations for a 3DOF full track chassis model
% https://www.overleaf.com/project/5e73ad44ea180c00018cb332 (Section 2.3.2)
%
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
%
% Author(s): 
% Tristan Pham       (atlpham@ucdavis.edu)        [Oct 2020 - ???     ]
% 
% Last Updated: 17-Feburary-2022

%% Test Case
if nargin == 0
    TFx = zeros(1,4);
    TFy = [250, 250, 0, 0];
    TMz = zeros(1,4);
    
    AFx = -100;
    AFy = 0;
    AMz = 0;
    
    Wheelbase  = 1.575; 
    TrackWidth = 1.22*ones(1,2);
    Steer      = [18, 15, 1, -1];
    
    Mass       = 275;   
    YawInertia = 130; 
    CoG        = [(0.47-0.5)*Wheelbase, 0, 0.25];
    
    LongVel = 15;
    LatVel  = 0;
    YawVel  = 0.5;
    
    [LongAcc, LatAcc, YawAcc, LongAccTot, LatAccTot] = ...
        ForceBasedRollCenter6DOFAccelerations( ...
            TFx, TFy, TMz, AFx, AFy, AMz, ...
            Wheelbase, TrackWidth, Steer, ... 
            Mass, YawInertia, CoG, ... 
            LongVel, LatVel, YawVel ) %#ok<NOPRT>
    
    return;
end

%% Computations

%%% Tire Positions (n,4,3 numeric)
TirePos = cat( 3, Wheelbase/2 + [-1, -1, 1, 1].*CoG(:,1), ...
                  [TrackWidth(:,1).*[1 -1]/2, TrackWidth(:,2).*[1 -1]/2], ...
                  zeros( size(Steer) ) );

%%% Tire Loads (n,4,3 numeric)
TireLoad = cat( 3, TFx, TFy, zeros( size(Steer) ) );

%%% Tire Yaw Moment (n,1 numeric)
TireMoment = cross(TirePos, TireLoad, 3);
TireMoment = sum( TireMoment(:,:,3) ) + sum(TMz,2);

%%% Moments
%Correct a or b
RollAcc = (sum((TFx.*sind(Steer) + TFy.*cosd(Steer)).*Front_View_Suspension_Arm...
    - (TFx.*cosd(Steer)-TFy.*sind(Steer)).*tand(Side_View_Jacking_Angle).*tw/2)...
    - AFy .* (Sprung_COG_Height + z) - AMz - (front_roll_stiffness + rear_roll_stiffness).*roll...
    - (roll_dampening_coeff_front + roll_dampening_coeff_rear).*roll_speed) / Ixx;

%Correct a or b
PitchAcc = (sum(-(TFx.*cosd(Steer) - TFy.*sind(Steer)).*Side_View_Suspension_Arm...
    + (TFx.*sind(Steer) + TFy.*cosd(Steer)).*tand(Front_View_Jacking_Angle).* (a or b) )...
    + AFx .* (Sprung_COG_Height + z) - AFz.*(a-L/2) - AMy ...
    - pitch_stiffness.*pitch - pitch_dampening.*pitch_velocity) ./ Iyy;

%Positive or negative trackwidth?
YawAcc = (sum((TFx.*cosd(Steer) - TFy.*sind(Steer)).*(-)tw/2 ...
    + (TFx.*sind(Steer) + TFy.*cosd(Steer)).*(a or -b) +TMz)...
    - AFy .* (L/2 - a) - AMz) ./ Izz; 

%%% Chassis Accelerations
LongAcc = ((sum( TFx.*cosd(Steer) - TFy.*sind(Steer), 2 ) - AFx) ./ Mass) + LatVel.*YawVel;

LatAcc  = ((sum( TFy.*cosd(Steer) + TFx.*sind(Steer), 2 ) - AFy) ./ Mass) - LongVel.*YawVel;

VertAcc = (sum(kr + br.*j + (TFx.*cosd(Steer) - TFy.*sind(Steer).*tand(Side_View_Jacking_Angle)...
    + (TFy.*cosd(Steer) + TFx.*sind(Steer)).*tan(Front_View_Jacking_Angle) - AFz) ./ Mass;

LongAccTot = LongAcc - LatVel .*YawVel;
LatAccTot  = LatAcc  + LongVel.*YawVel;

%%% Matricies (FIXME)
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