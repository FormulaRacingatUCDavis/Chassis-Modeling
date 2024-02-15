
function [ddzs, ddzu, ddzt] = quarterCarModel(zs, zu, zt, dzs, dzu)
%basics - calculating the accelerations of sprung mass, unsprung mass, and
%tire of 1/4 of a car. 

    % parameters (arbitrary values)
    ms = 50; % sprung mass in kg (mass of the car body)
    mu = 10;  % unsprung mass in kg (mass of the wheel assembly)
    mt = 4;  % tire mass in kg
    ks = 50000; % spring constant of suspension in N/m
    kt = 127000; % spring constant of tire in N/m
    d = 0.5;   % damping coefficient in Ns/m

    % dynamics equations
    % sprung mass (ddzs)
    ddzs = (ks * ((zu - zs)+0.5) - d * (dzu - dzs)) / ms;
    %note z_ is merely referring to the vertical positions of the
    %components

    % unsprung mass (ddzu)
    ddzu = (ks * ((zs - zu)-0.5) + (kt * ((zt - zu)+0.5)) - d * (dzs - dzu)) / mu;

    % tire (ddzt)
    ddzt = kt * ((zu - zt)-0.5) / mt;

end

%in live script, 
% test inputs
%zs_test = 0.1;
%zu_test = 0.05;
%zt_test = 0;
%[ddzs, ddzu, ddzt] = quarterCarModel(zs_test, zu_test, zt_test)