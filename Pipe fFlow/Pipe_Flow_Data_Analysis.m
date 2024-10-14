%{
MECH&AE 157W - Basic Mechanical and Aerospace Engineering Laboratory with Writing
Instructor: Prof. Yongie Hu
Teaching Assistants: Zihao Qin, Qiyu Xing, Benjamin Heronimus
Fall 2024 – UCLA
Students: Alex Lie, Andrew Tan, Anli Liu, Umer Badae, Ian Lee
Lab Group: 6
Lab: Pipe Flow
%}

%Beginning of code
clear
close
clc

%Importing data
data=readtable("Pipe_Flow_Data.xlsx");
Time = data.Time; %Units: s
PipeNumber = data.PipeNumber; %Unitless
Temperature = data.Temperature; %Units: Kelvin
Frequency = data.Frequency; %Units: Hz
PressureDrop = data.PressureDrop; %Units: Pa

%Pipe Inner Diamaters (in order of 1, 2, 3, 4)
PipeDiameterArray = [8.15; 11.8; 17.5; 9.93] ./ 1000; %Units: m

%Pipe Lengths (in order of 1, 2, 3, 4)
PipeLengthArray = [0.781; 1.143; 1.98; 0.337]; %Units: m

% create arrays of length 40 with associated pipe diameter
PipeDiameter = zeros(40,1);
PipeLength = zeros(40,1);
for a=0:3
    PipeDiameter(10*a+1:10*(a+1)) = PipeDiameterArray(a+1);
    PipeLength(10*a+1:10*(a+1)) = PipeLengthArray(a+1);
end
PipeCrossSectionalArea= pi*PipeDiameter.^2/4; %Units: m^2

k = 5244;
Velocity = ((Frequency ./ k) .* 0.00378541) ./ ((pi .* PipeDiameter.^2) ./ 4); % m/s

Density = zeros(40,1);
DynViscosity = zeros(40,1);
for a=1:40
    Density(a) = waterDensity(Temperature(a));
    DynViscosity(a) = waterDynViscosity(Temperature(a));
end

Re = (Density .* Velocity .* PipeDiameter) ./ DynViscosity;

f_theoretical = (0.790 * log(Re) - 1.64).^(-2); % theoretical empirical calculation for smooth pipe
h = 0.305 / 1000;
P = 3.08 / 1000;
epsilon_s = h * exp(3.4 - 0.42 * (P / h)^0.46);
f_theoretical(31:40) = (1.74 + 2 * log10(PipeDiameter(31:40) / (2 * epsilon_s))).^(-2); % theoretical empirical calculation for rough pipe

f_experimental = PressureDrop ./ (((PipeLength ./ PipeDiameter)) .* (0.5 * Density .* (Velocity).^2));

%% poly (1st degree) fitting

p1 = polyfit(Re(1:9), f_experimental(1:9), 1);
x1 = linspace(Re(1), Re(10), 10);
y1 = polyval(p1, x1);

p2 = polyfit(Re(11:19), f_experimental(11:19), 1);
x2 = linspace(Re(11), Re(20), 10);
y2 = polyval(p2, x2);

p3 = polyfit(Re(21:29), f_experimental(21:29), 1);
x3 = linspace(Re(21), Re(30), 10);
y3 = polyval(p3, x3);

p4 = polyfit(Re(31:39), f_experimental(31:39), 1);
x4 = linspace(Re(31), Re(40), 10);
y4 = polyval(p4, x4);

%%
% plot stuff
figure(1);
hold on
plot(Re(1:10), f_experimental(1:10), '-o');
plot(Re(1:10), f_theoretical(1:10));
plot(x1,y1);
hold off
xlabel('Re');
ylabel('Friction Factor');
title('Moody Plot for Small Smooth Pipe');
legend('Experimental', 'Theoretical', 'Fitted Line');

figure(2);
hold on
plot(Re(11:20), f_experimental(11:20), '-o');
plot(Re(11:20), f_theoretical(11:20));
plot(x2,y2);
hold off
xlabel('Re');
ylabel('Friction Factor');
title('Moody Plot for Medium Smooth Pipe');
legend('Experimental', 'Theoretical', 'Fitted Line');

figure(3);
hold on
plot(Re(21:30), f_experimental(21:30), '-o');
plot(Re(21:30), f_theoretical(21:30));
plot(x3,y3);
hold off
xlabel('Re');
ylabel('Friction Factor');
title('Moody Plot for Large Smooth Pipe');
legend('Experimental', 'Theoretical', 'Fitted Line');

figure(4);
hold on
plot(Re(31:40), f_experimental(31:40), '-o');
plot(Re(31:40), f_theoretical(31:40));
plot(x4,y4);
hold off
xlabel('Re');
ylabel('Friction Factor');
title('Moody Plot for Rough Pipe');
legend('Experimental', 'Theoretical', 'Fitted Line');
%% Formulas
%% Table approximations
% 1st order approximation for water density (kg/m^3)
function ans = waterDensity(temperatureK)
    m = -0.00012;
    b = 1;
    ans = (b + m * (temperatureK - 273.15)) * 1000;
end

% 1st order approximation for water dynamic viscosity (kg / m*s)
function ans = waterDynViscosity(temperatureK)
    m = -2.05e-5;
    b = 0.891e-3;
    ans = b + m * (temperatureK - 273.15 - 25);
end

%% calculate friction factor
% experimental
function f = experimentalFriction(Pressure_drop, Diameter, Length, density, flow_velocity)
    f = (Pressure_drop * Diameter / Length)/(0.5 * density * flow_velocity^2);
end

% theoretical laminar
function f = theoreticalLaminarFrictionSmooth(Re)
    f = 64/Re;
end

% theoretical turbulent smooth - Petukhov formula - valid for 1E4 < Re < 5E6
function f = theoreticalTurbulentFrictionSmooth(Re_D)
    f = (0.790 * ln(Re_D) - 1.64)^-2;
end

% theoretical turbulent rough - Nikuradse's formula
function f = theoreticalTurbulentFrictionRough()
    % There's only one rough pipe, so the roughness is directly calculated inside this function
    D = 9.93E-3;
    h = 0.305E-3;
    P = 3.08E-3;
    roughness = h * exp(3.4 - 0.42 * (P/h)^0.46);
    f = (1.74 + 2.0*log10(0.5*D/roughness))^-2;
end

%% Error Analysis
% standard error
function Sy = StandardError(Y_actual, Y_fitted)
    N = length(Y_actual);
    sum = 0;
    for i = 1:1:N
        sum = sum + (Y_actual - Y_fitted)^2;
    end
    Sy = sqrt(sum / (N - 2));
end

% Precision of a measurement
function Py = Precision_measurement(Y_actual, Y_fitted)
    % spits out a vector with each of the precision uncertainty for each measurement
    N = length(Y_actual);
    x_avg = sum(Y_actual)/N;
    sum_square = 0;
    for i = 1:1:N
        sum_square = sum_square + (Y_actual(i))^2;
    end
    Sxx = sum_square - 1/N * sum(Y_actual)^2;
    Py = zeros(N,1);
    for i = 1:1:N
        Py(i) = 2 * sqrt(StandardError(Y_actual, Y_fitted)^2 * (1 + 1/N + (Y_actual(i)-x_avg)^2/Sxx));
    end
end

% Precision of a curve
function Py = Precision_curve(Y_actual, Y_fitted)
    % spits out a vector with each of the precision uncertainty for each measurement
    N = length(Y_actual);
    x_avg = sum(Y_actual)/N;
    sum_square = 0;
    for i = 1:1:N
        sum_square = sum_square + (Y_actual(i))^2;
    end
    Sxx = sum_square - 1/N * sum(Y_actual)^2;
    Py = zeros(N,1);
    for i = 1:1:N
        Py(i) = 2 * sqrt(StandardError(Y_actual, Y_fitted)^2 * (1/N + (Y_actual(i)-x_avg)^2/Sxx));
    end
end