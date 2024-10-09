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
PipeDiameter = [8.15; 11.8; 17.5; 9.93] ./ 1000; %Units: m
PipeCrossSectionalArea= pi*PipeDiameter.^2/4; %Units: mm^2

%Pipe Lengths (in order of 1, 2, 3, 4)
PipeLength = [0.781; 1.143; 1.98; 0.337]; %Units: m

k = 5244;
Velocity = Frequency ./ k;

Density = zeros(40,1);
DynViscosity = zeros(40,1);
for a=1:40
    Density(a) = waterDensity(Temperature(a));
    DynViscosity(a) = waterKinViscosity(Temperature(a));
end

Re1 = (Density(1:10) .* Velocity(1:10) .* PipeDiameter(1)) ./ DynViscosity(1:10);
Re2 = (Density(11:20) .* Velocity(11:20) .* PipeDiameter(2)) ./ DynViscosity(11:20);
Re3 = (Density(21:30) .* Velocity(21:30) .* PipeDiameter(3)) ./ DynViscosity(21:30);
Re4 = (Density(31:40) .* Velocity(31:40) .* PipeDiameter(4)) ./ DynViscosity(31:40);

f1_theoretical = (0.790 * log(Re1) - 1.64).^(-2);
f2_theoretical = (0.790 * log(Re2) - 1.64).^(-2);
f3_theoretical = (0.790 * log(Re3) - 1.64).^(-2);
f4_theoretical = (0.790 * log(Re4) - 1.64).^(-2);

f1 = (PressureDrop(1:10) .* (Velocity(1:10) ./ PipeLength(1))) ./ (0.5 * Density(1:10) .* (Velocity(1:10)).^2);
f2 = (PressureDrop(11:20) .* (Velocity(11:20) ./ PipeLength(2))) ./ (0.5 * Density(11:20) .* (Velocity(11:20)).^2);
f3 = (PressureDrop(21:30) .* (Velocity(21:30) ./ PipeLength(3))) ./ (0.5 * Density(21:30) .* (Velocity(21:30)).^2);
f4 = (PressureDrop(31:40) .* (Velocity(31:40) ./ PipeLength(4))) ./ (0.5 * Density(31:40) .* (Velocity(31:40)).^2);

% plot stuff

figure(1);
hold on
%plot(Re1, f1, '-o');
plot(Re1, f1_theoretical);
hold off
xlabel('Re');
ylabel('Friction Factor');
title('Moody Plot for Small Smooth Pipe');

figure(2);
plot(Re2, f2, '-o');
xlabel('Re');
ylabel('Friction Factor');
title('Moody Plot for Medium Smooth Pipe');

figure(3);
plot(Re3, f3, '-o');
xlabel('Re');
ylabel('Friction Factor');
title('Moody Plot for Large Smooth Pipe');

figure(4);
plot(Re4, f4, '-o');
xlabel('Re');
ylabel('Friction Factor');
title('Moody Plot for Rough Pipe');

% 1st order approximation for water density (kg/m^3)
function ans = waterDensity(temperatureK)
    m = -0.00012;
    b = 1;
    ans = b + m * (temperatureK - 273.15);
end

% 1st order approximation for water dynamic viscosity (Pa*s)
function ans = waterKinViscosity(temperatureK)
    m = -0.0253e-6;
    b = 1.004e-6;
    ans = b + m * (temperatureK - 273.15 - 20);
end