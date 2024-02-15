function [ddpsi, ddx, ddy] = bicycleModelDynamics(dpsi, dx, dy, delta)
    % Parameters (arbitrary values)
    L = 1.525; % distance from CG to front wheel
    I = 130 % vehicle moment of inertia
    m = 275; % mass of the vehicle
    
    % Forces
    Fry = 1000; 
    Ffx = 1000;
    Ffy = 2000; 
    Frx = 2000;

    % Dynamics equations
    Lf = L*(1-PercentFront);
    Lr = L*(PercentFront);
    ddpsi = (Lr * Fry - Lf * (Ffx * sin(delta) + Ffy * sin(delta))) / I;
    ddx = (Frx + Ffx * cos(delta) - Ffy * sin(delta) - m * dpsi * dy) / m;
    ddy = (Fry + Ffx * sin(delta) + Ffy * sin(delta) - m * dpsi * dx) / m;
    
    % Output the second derivatives 
end

%current_dpsi = 0.1 % yaw rate in radians per sec
%current_dx = 10 % longitudinal velocity m/s
%current_dy = 10 % lateral velocity m/s
%current_delta = 10 % steering angle in radians

%[ddpsi, ddx, ddy] = bicycleModelDynamics(current_dpsi, current_dx, current_dy, current_delta);