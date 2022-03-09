function [TireLateralForce] = tireModel (gamma, alpha, Fz)
    a = [1.3 -7*10^-2 1.1 1.18 7.8 0 -0.2 2.4*10^-2 2.53*10^-2 0 0 2.57*10^-3 0 0];
    C = a(1);
    Sh = a(9)*gamma + a(10)*Fz + a(11);              %horizontal shift
    Sv = a(12)*Fz*gamma + a(13)*Fz + a(14);          %vertical shift
    x = alpha + Sh;                                  %slip angle
    D = (a(2)*Fz + a(3))*Fz;
    BCD = a(4)*(sin(2*atan(Fz/a(5))))*(1 - a(6)*gamma);
    B = BCD/(C*D);
    E = a(7)*Fz + a(8);
    TireLateralForce = D*sin(C*atan(B*x - E(B*x - atan(B*x)))) + Sv;
end
