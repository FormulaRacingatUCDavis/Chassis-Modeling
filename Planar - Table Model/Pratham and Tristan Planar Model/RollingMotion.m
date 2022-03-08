function [RollAngle, FrontLoadShift, RearLoadShift] = RollingMotion(FrontRollStiffness,...
    RearRollStiffness, y_ddot, FrontDampeningCoeff, RearDampeningCoeff, )

%% RollingMotion - Roll Evaluation
% This model evaluates rolling motion of a planar vehicle model 
% 
% Parameters:
%   SprungMass         - (n,1 numeric) Sprung Mass                           {m_s} [kg]
%   hs                 - (n,1 numeric) Height from roll center to COG        {hs} [m]
%   g                  - (n,1 numeric) Acceleration due to gravity           {g} [m/s^2]
%   tw                 - (n,1 numeric) Track width                           {tw} [m]
%   hf                 - (n,1 numeric) Front height from road to roll center {hf} []
%   hr                 - (n,1 numeric) Rear height from road to roll center  {hr} []
%
% Inputs:
%   FrontRollStiffness  - (n,1 numeric) Front roll stiffness        {K_(phi_f)} []
%   RearRollStiffness   - (n,1 numeric) Rear roll stiffness         {K_(phi_r)} []
%   RollAngle           - (n,1 numeric) Roll angle of vehicle body  {phi} [deg]
%   y_ddot              - (n,1 numeric) Lateral Acc                 {y_ddot} [m/s^2]   
%   FrontDampeningCoeff - (n,1 numeric) Front Dampening Coefficient {C_phif} []
%   RearDampeningCoeff  - (n,1 numeric) Front Dampening Coefficient {C_phir} []
%   RollVelocity        - (n,1 numeric) Front Dampening Coefficient {phi_dot} []
%
% Outputs:
%   FrontLoadShift - (n,1 numeric) Front Load transfer {deltaF_zf} [N]
%   RearLoadShift  - (n,1 numeric) Rear Load transfer  {deltaF_zr} [N]
%   RollAngle      - (n,1 numeric) Rear Load transfer  {phi} [degrees]
%
% Notes:
% 
%
% Author(s): 
% Tristan Pham       (atlpham@ucdavis.edu       ) [Jan 2021 -         ]
% 
% Last Updated: 5-March-2022

%% Test Case

%% Computation

%%%Parameters
SprungMass = ;
g = 9.8;
hs =;
hf = ;
hr = ;
tw = ;

%%%Equation of rolling motion (FIXME HOW DO I SOLVE FOR ROLL ANGLE)
RollAngle = (y_ddot*SprungMass*hs + ...
    SprungMass*g*hs*sin(RollAngle)) / (FrontRollStiffness+RearRollStiffness) ;

%%%Load shift equations
FrontLoadShift = (FrontRollStiffness*RollAngle + FrontDampeningCoeff*RollVelocity +...
    (Fy1 + Fy2)*hf) * 1/tw;
RearLoadShift = (RearRollStiffness*RollAngle + RearDampeningCoeff*RollVelocity +...
    (Fy3 + Fy4)*hr) * 1/tw;
    
