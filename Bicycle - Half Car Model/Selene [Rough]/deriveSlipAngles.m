function [alphaF, alphaR] = deriveSlipAngles(delta, omega,dx)    
    % Parameters:
    L = 1.525; 
    PercentFront = 0.6;

    % Inputs:
    % delta  - Steering angle (radians)
    % omega  - Yaw rate (radians per second)
    % dx     - Longitudinal velocity (m/s)
    
    % Outputs:
    % alphaF - Front slip angle (radians)
    % alphaR - Rear slip angle (radians)

    Lf = L*(1-PercentFront);
    Lr = L*(PercentFront);

    % Derive the lateral velocity dy
    % Assuming small slip angles and a steady-state cornering condition
    dy = omega * (Lf + Lr); % Since v = R * omega and R = a + b for small slip angles

    % Using small angle approximations for the arctan function
    alphaF = delta - (dy + Lf * omega) / dx;
    alphaR = (dy - Lr * omega) / dx;
end
