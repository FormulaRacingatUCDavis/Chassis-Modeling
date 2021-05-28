function [RearSteeringAngle, LongAccTot, LatAccTot, YawAcc] = ...
    Acceleration_Function(SteeringWheelAngle,LongForce,LatForce,...
    Mass,WheelBase,PercentFront,TrackWidth)
%% Planar Model - Planar Model Acceleration Calculations
% 
% Inputs:
%  SteeringWheelAngle - (n,1 numeric) SteeringWheelAngle {delta} [degrees]
%  LongForce          - (n,1 numeric) LongitudinalForce  {Fx} [N-m]
%  LatForce           - (n,1 numeric) Lateral Force      {Fy} [N-m]
%  Mass               - (n,1 numeric) Mass               {m} [kg]
%  WheelBase          - (n,1 numeric) Wheel Base         {L} [m]
%  PercentFront       - (n,1 numeric) Percent Front      {pf} [%]
%  TrackWidth         - (n,1 numeric) TrackWidth         {tw} [m]
%
% Outputs:
%  RearSteeringAngle - (n,1 numeric) Rear Steering Angle             {delta} [degrees]
%  LongAccTot        - (n,1 numeric) Total Longitudinal Acceleration {a_x} [m/s^2]
%  LatAccTot         - (n,1 numeric) Total Lateral Acceleration      {a_y} [m/s^2]
%  YawAcc            - (n,1 numeric) Yaw Acceleration                {psi_ddot} [rad/s^2]
%
% Notes:
 
% 
%
% Author(s): 
% Tristan Pham (atlpham@ucdavis.edu) [Sep 2020 - Jun 2021] 

% Last Updated: 24-May-2021


%% Test Cases
if nargin == 0
    %%% Test Inputs
     SteeringWheelAngle = ; 
     LongForce= ; 
     LatForce = ; 
     Mass = ;
     WheelBase = ;
     PercentFront = ;
     TrackWidth = ;
    
    fprintf('Executing Planar Model() Test Cases: \n');
    
    
    [RearSteeringAngle, LongAccTot, LatAccTot, YawAcc] = BrakeModel(SteeringWheelAngle,...
        LongForce,LatForce,Mass,WheelBase,PercentFront,TrackWidth);
    
 % for i = 1:numel()
      %  fprintf('   Steady State Instance %i: \n', i);
      % fprintf('      tau_i = %5.2f [N-m] \n', InputTorque(i));
   % end
    
    
     %return;   
end
    
%% Computation
% Parameters 
TrackWidth = TrackWidth/39.3700787;
WheelBase = WheelBase/39.3700787;  
Mass = Mass/2.20462262;     %Conversions
a = WheelBase * (1 - PercentFront);     %Because the ceneter or mass is more to there rear so a is longer
b = WheelBase * PercentFront;
YawInertia = Mass * (1.3)^2;  %Cacls

% Steer Angle Calculation
RearSteeringAngle = 0.25 * SteeringWheelAngle;  %Steering angle ratio to the rear steergin angle 

% Logitudinal
LongAccTot = ( LongForce(1).*cosd(RearSteeringAngle) + LongForce(2).*cosd(RearSteeringAngle) - ... 
        LatForce(1).*sind(RearSteeringAngle) - LatForce(2).*sind(RearSteeringAngle) + LongForce(3) + LongForce(4) ) / Mass; 

% Lateral
LatAccTot = ( LatForce(1).*cosd(RearSteeringAngle) + LongForce(1).*sind(RearSteeringAngle) + ...
        LatForce(2).*cosd(RearSteeringAngle) + LongForce(2).*sind(RearSteeringAngle) + LatForce(3) + LatForce(4) ) / Mass;

% Yaw
YawAcc = ( a.*( LatForce(1) + LatForce(2) ) - b.*( LatForce(3) + LatForce(4) ) + ...
            (TrackWidth/2).*( LongForce(2) + LongForce(4) ) - (TrackWidth/2).*( LongForce(1) + LongForce(3) ) ) / YawInertia;