function [TireLateralForce] = PlanarTireModel (gamma, alpha, Fz)
    %% Tire Force - Evaluation of Lateral Tire Force
% This model evaluates lateral tire force of a planar vehicle model using 
% Pacejka Magic Formula.
% 
% Parameters:
%   Magic Formula Constant      - (n,14 numeric) Sprung Mass                 {m_s} [kg]
%
% Inputs:
%   gamma                       - (n,1 numeric)  Front roll stiffness        {gamma} [deg]
%   alpha                       - (n,1 numeric)  Slip Angle         {alpha)} [deg]
%   Fz                          - (n,1 numeric)  Roll angle of vehicle body  {Fz} [Newtons]
%
% Outputs:
%   TireLateralForce            - (n,1 numeric) Tire Lateral Force           {Y} [N]
%
% Notes:
% 
%
% Author(s): 
% Tristan Pham       (atlpham@ucdavis.edu       ) [Jan 2021 -         ]
% 
% Last Updated: 12-March-2022

%% Test Case

%% Computation

%%%Parameters
    a = [1.3 -7*10^-2 1.1 1.18 7.8 0 -0.2 2.4*10^-2 2.53*10^-2 0 0 2.57*10^-3 0 0];
    
%%% Equations
    C = a(1);                                        % Shape Factor
    Sh = a(9)*gamma + a(10)*Fz + a(11);              % Horizontal Shift
    Sv = a(12)*Fz*gamma + a(13)*Fz + a(14);          % Vertical Shift
    x = alpha + Sh;                                  % Slip Angle
    D = (a(2)*Fz + a(3))*Fz;                         % Peak Value
    BCD = a(4)*(sin(2*atan(Fz/a(5))))*(1 - a(6)*gamma);
    B = BCD/(C*D);                                   % Stiffness Factor
    E = a(7)*Fz + a(8);                              % Curvature Factor
    TireLateralForce = D*sin(C*atan(B*x - E(B*x - atan(B*x)))) + Sv;
end
