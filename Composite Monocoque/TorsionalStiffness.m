%% Preambles
clc
clear
close all

%% User Input
XDistance_mm = 100:100:2000;  % X distance from rear axle.

YDistance_mm = [...  % Y cross-sectional region width to midplane.
    261 265.7056 298.063 330.4893 362.9155 368.95 368.5464 364.5259 ...
    356.174 343.4434 326.2607 304.5244 278.1014 248.7583 235.6 235.6...
    235.5976 229.6597 215.201 194.8652];

ZDistance_mm = [...  % Z displacement due to loading.
    0.029999 0.074716 0.15963 0.24602 0.33503 0.3813 0.4157 0.44194 ...
    0.45739 0.46083 0.45096 0.42647 0.38493 0.33275 0.34221 0.3951  ...
    0.43956 0.46465 0.4643 0.44709];

TrackWidth_mm = 1153.3846;

Load_N = 2000;

%% Data Processing
XDistance_m = XDistance_mm ./ 1000;
YDistance_m = YDistance_mm ./ 1000;
ZDistance_m = ZDistance_mm ./ 1000;
TrackWidth_m = TrackWidth_mm ./ 1000;
TwistAngle_deg = atand(ZDistance_m ./ YDistance_m);

[cubx, cuby] = CUBIC2(XDistance_m, TwistAngle_deg);  % Natural cubic spline

TorsionalStiffness_Nm_deg = (Load_N*TrackWidth_m) ./ cuby;

%% Plotting
figure;
hold on;
plot(cubx, cuby, '.b');
scatter(XDistance_m, TwistAngle_deg, '^r');
title("Twist Angle Along Length of Monocoque");
xlabel('Distance from Rear Axle (m)');
ylabel('Twist Angle (deg)');
legend('Natural Cubic Spline', 'Raw Data Points');

figure;
plot(cubx, TorsionalStiffness_Nm_deg);
title({"Torsional Stiffness Along Length of Monocoque", ...
    "Minimum TS: " + min(TorsionalStiffness_Nm_deg)});
xlabel('Distance from Rear Axle (m)');
ylabel('Torsional Stiffness (Nm/deg)');

%% Clear Extra Variables
clear cubx cuby delta dTwistAngle_deg_m i Load_N ...
    TorsionalStiffness_Nm_deg TrackWidth_m TrackWidth_mm TwistAngle_deg ...
    XDistance_m XDistance_mm YDistance_m YDistance_mm ZDistance_m ...
    ZDistance_mm

%% Supporting Functions
function [cubx, cuby] = CUBIC2(x, y)
% Natural cubic spline function.

N = length(x);
delta_x = x(2) - x(1);  % Assume uniform.

% Generate tridiagonal matrix.
a = zeros(1, N);
a(1) = 0;
a(2:end-1) = 1;
a(end) = 0;
b = zeros(1, N);
b(1) = 1;
b(2:end-1) = 4;
b(end) = 1;
c = zeros(1, N);
c(1) = 0;
c(2:end-1) = 1;
c(end) = 0;
d = zeros(1, N);
d(1) = 0;
for i = 1:length(d)-2
    d(i+1) = (6 * (y(i+2) - 2*y(i+1) + y(i))) / (delta_x^2);
end
d(end) = 0;

m = THOMAS3(a, b, c, d);  % Solve tridiagonal matrix.
a = zeros(1, N-1);
b = (1/2) .* m;
c = zeros(1, N-1);
d = y(1:N-1);

for i = 2:length(m)
    delta_y = y(i) - y(i-1);
    a(i-1) = (m(i) - m(i-1)) / (6*delta_x);
    c(i-1) = (delta_y/delta_x) - (m(i-1)*delta_x/3) - (m(i)*delta_x/6);
end

cubx = [];
cuby = [];

% Using generated coefficients to interpolate values.
for i = 1:length(x)-1
    N = 10;
    delta = (x(i+1)-x(i))/N;
    tempx = linspace(x(i), x(i+1)-delta, N);  % 10 values between each node.
    tempy = [];
    for k = 1:length(tempx)
        tempy = [tempy, a(i)*(tempx(k)-x(i))^3 + ...
            b(i)*(tempx(k)-x(i))^2 + c(i)*(tempx(k)-x(i)) + d(i)];
    end
    cubx = [cubx, tempx];
    cuby = [cuby, tempy];
end

end

function [x] = THOMAS3(a, b, c, d)
% Tridiagonal matrix solver.

n = length(b);

for i = 1:n-1
    b(i+1) = b(i+1) - (c(i) * a(i+1) / b(i));
    d(i+1) = d(i+1) - (d(i) * a(i+1) / b(i));
    a(i+1) = 0;
    c(i) = c(i)/b(i);
    d(i) = d(i)/b(i);
    b(i) = 1;
end

d(n) = d(n)/b(n);
b(n) = 1;
x = zeros(1,n);
x(n) = d(n);
for i = n-1:-1:1
    x(i) = d(i) - (c(i) * x(i+1));
end

end
