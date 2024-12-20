%{
MECH&AE 157W - Basic Mechanical and Aerospace Engineering Laboratory with Writing
Instructor: Prof. Yongie Hu
Teaching Assistants: Zihao Qin, Qiyu Xing, Benjamin Heronimus
Fall 2024 – UCLA
Students: Alex Lie, Andrew Tan, Anli Liu, Umer Badar, Ian Lee
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
f_theoretical(31:39) = (1.74 + 2 * log10(PipeDiameter(31:39) / (2 * epsilon_s))).^(-2); % theoretical empirical calculation for rough pipe

f_experimental = PressureDrop ./ (((PipeLength ./ PipeDiameter)) .* (0.5 * Density .* (Velocity).^2));

%% poly (1st degree) fitting

p1 = polyfit(Re(1:10), log(f_experimental(1:10)), 1);% Log fit for the error
y1 = polyval(p1, Re(1:10));

p2 = polyfit(Re(11:19), log(f_experimental(11:19)), 1);% Log fit for the error
y2 = polyval(p2, Re(11:19));

p3 = polyfit(Re(21:29), log(f_experimental(21:29)), 1);% Log fit for the error
y3 = polyval(p3, Re(21:29));

p4 = polyfit(Re(31:39), log(f_experimental(31:39)), 1);% Log fit for the error
y4 = polyval(p4, Re(31:39));

%% Error Analysis 

% Pipe 1
% Calculate the standard error 
Sy1 = sqrt(1 / (length(Re(1:10)) - 2) * sum((f_experimental(1:10) - exp(y1)).^2)); % Standard error
xAvg1 = mean(Re(1:10)); % Average of Reynolds numbers
Sxx1 = sum((Re(1:10)).^2) - (1 / length(Re(1:10))) * (sum(Re(1:10)))^2; % Sample variance

% Calculate Py and Py_hat 
Py1 = 2 * sqrt(Sy1.^2 * (1 + 1 / length(Re(1:10)) + (Re(1:10) - xAvg1).^2 ./ Sxx1)); % Precision uncertainty of measurement
Py_hat1 = 2 * sqrt(Sy1.^2 * (1 / length(Re(1:10)) + (Re(1:10) - xAvg1).^2 ./ Sxx1)); % Precision uncertainty of curve fit


% Pipe 2
% Calculate the standard error
Sy2 = sqrt(1 / (length(Re(11:19)) - 2) * sum((f_experimental(11:19) - exp(y2)).^2)); % Standard error
xAvg2 = mean(Re(11:19)); % Average of Reynolds numbers
Sxx2 = sum((Re(11:19)).^2) - (1 / length(Re(11:19))) * (sum(Re(11:19)))^2; % Sample variance
% Calculate Py and Py_hat 
Py2 = 2 * sqrt(Sy2.^2 * (1 + 1 / length(Re(11:19)) + (Re(11:19) - xAvg2).^2 ./ Sxx2)); % Precision uncertainty of measurement
Py_hat2 = 2 * sqrt(Sy2.^2 * (1 / length(Re(11:19)) + (Re(11:19) - xAvg2).^2 ./ Sxx2)); % Precision uncertainty of curve fit


%Pipe 3
% Calculate the standard error
Sy3 = sqrt(1 / (length(Re(21:29)) - 2) * sum((f_experimental(21:29) - exp(y3)).^2)); % Standard error
xAvg3 = mean(Re(21:29)); % Average of Reynolds numbers
Sxx3 = sum((Re(21:29)).^2) - (1 / length(Re(21:29))) * (sum(Re(21:29)))^2; % Sample variance
% Calculate Py and Py_hat 
Py3 = 2 * sqrt(Sy3.^2 * (1 + 1 / length(Re(21:29)) + (Re(21:29) - xAvg3).^2 ./ Sxx3)); % Precision uncertainty of measurement
Py_hat3 = 2 * sqrt(Sy3.^2 * (1 / length(Re(21:29)) + (Re(21:29) - xAvg3).^2 ./ Sxx3)); % Precision uncertainty of curve fit


%Pipe 4
% Calculate the standard error
Sy4 = sqrt(1 / (length(Re(31:39)) - 2) * sum((f_experimental(31:39) - exp(y4)).^2)); % Standard error
xAvg4 = mean(Re(31:39)); % Average of Reynolds numbers
Sxx4 = sum((Re(31:39)).^2) - (1 / length(Re(31:39))) * (sum(Re(31:39)))^2; % Sample variance
% Calculate Py and Py_hat 
Py4 = 2 * sqrt(Sy4.^2 * (1 + 1 / length(Re(31:39)) + (Re(31:39) - xAvg4).^2 ./ Sxx4)); % Precision uncertainty of measurement
Py_hat4 = 2 * sqrt(Sy4.^2 * (1 / length(Re(31:39)) + (Re(31:39) - xAvg4).^2 ./ Sxx4)); % Precision uncertainty of curve fit














% plot stuff
                        %pipe 1
figure(1);
hold on
loglog(Re(1:10), f_experimental(1:10), 'x','Color', [0 1 0], 'lineWidth', 2); % experimental friction
loglog(Re(1:10), f_theoretical(1:10), 'o','Color', [1 0 0], 'lineWidth', 2); % theoretical friction
plot(Re(1:10),exp(y1),'Color', [0 0 1], 'lineWidth', 2); % Fitted Line

hold off
xlabel('log(Re)');
ylabel('log(fr) [Friction Factor]');
title('Moody Plot for Small Smooth Pipe');
legend('Experimental Data', 'Theoretical (Petukhov) Data', 'Fitted Line');
                  
figure(2);
hold on
errorbar(Re(1:10), f_experimental(1:10),Py1, 's','LineWidth',1.5); % error bar
plot(Re(1:10),exp(y1),'Color', [0 0 1], 'lineWidth', 2); % Fitted Line
plot(Re(1:10),exp(y1)+Py_hat1,'--','Color', [0 0 0], 'lineWidth', 2);
plot(Re(1:10),exp(y1)-Py_hat1,'--','Color', [0 0 0], 'lineWidth', 2);

hold off
xlabel('log(Re)');
ylabel('log(fr) [Friction Factor]');
title('Moody Plot for Small Smooth Pipe');
legend('Error Bar (Precision uncertainty of a measurement)', 'Fitted Line','Error Lines (Precision uncertainity of the curve fit)');
                        %pipe 2
figure(3);
hold on
loglog(Re(11:19), f_experimental(11:19), 'x','Color', [0 1 0], 'lineWidth', 2); % experimental friction
loglog(Re(11:19), f_theoretical(11:19), 'o','Color', [1 0 0], 'lineWidth', 2); % theoretical friction
plot(Re(11:19),exp(y2),'Color', [0 0 1], 'lineWidth', 2); % Fitted Line

hold off
xlabel('log(Re)');
ylabel('log(fr) [Friction Factor]');
title('Moody Plot for Medium Smooth Pipe');
legend('Experimental Data', 'Theoretical (Petukhov) Data', 'Fitted Line');

figure(4);
hold on
errorbar(Re(11:19), f_experimental(11:19),Py2, 's','LineWidth',1.5); % error bar
plot(Re(11:19),exp(y2),'Color', [0 0 1], 'lineWidth', 2); % Fitted Line
plot(Re(11:19),exp(y2)+Py_hat2,'--','Color', [0 0 0], 'lineWidth', 2);
plot(Re(11:19),exp(y2)-Py_hat2,'--','Color', [0 0 0], 'lineWidth', 2);

hold off
xlabel('log(Re)');
ylabel('log(fr) [Friction Factor]');
title('Moody Plot for Medium Smooth Pipe');
legend('Error Bar (Precision uncertainty of a measurement)', 'Fitted Line','Error Lines (Precision uncertainity of the curve fit)');

                         %pipe 3
figure(5);
hold on
loglog(Re(21:29), f_experimental(21:29), 'x','Color', [0 1 0], 'lineWidth', 2); % experimental friction
loglog(Re(21:29), f_theoretical(21:29), 'o','Color', [1 0 0], 'lineWidth', 2); % theoretical friction
plot(Re(21:29),exp(y3),'Color', [0 0 1], 'lineWidth', 2); % Fitted Line

hold off
xlabel('log(Re)');
ylabel('log(fr) [Friction Factor]');
title('Moody Plot for Large Smooth Pipe');
legend('Experimental Data', 'Theoretical (Petukhov) Data', 'Fitted Line');

figure(6);
hold on
errorbar(Re(21:29), f_experimental(21:29),Py3, 's','LineWidth',1.5); % error bar
plot(Re(21:29),exp(y3),'Color', [0 0 1], 'lineWidth', 2); % Fitted Line
plot(Re(21:29),exp(y3)+Py_hat3,'--','Color', [0 0 0], 'lineWidth', 2);
plot(Re(21:29),exp(y3)-Py_hat3,'--','Color', [0 0 0], 'lineWidth', 2);

hold off
xlabel('log(Re)');
ylabel('log(fr) [Friction Factor]');
title('Moody Plot for Large Smooth Pipe');
legend('Error Bar (Precision uncertainty of a measurement)', 'Fitted Line','Error Lines (Precision uncertainity of the curve fit)');


                         %pipe 4
figure(7);
hold on
loglog(Re(31:39), f_experimental(31:39), 'x','Color', [0 1 0], 'lineWidth', 2); % experimental friction
loglog(Re(31:39), f_theoretical(31:39), 'o','Color', [1 0 0], 'lineWidth', 2); % theoretical friction
plot(Re(31:39),exp(y4),'Color', [0 0 1], 'lineWidth', 2); % Fitted Line

hold off
xlabel('log(Re)');
ylabel('log(fr) [Friction Factor]');
title('Moody Plot for Square Rough Pipe');
legend('Experimental Data', 'Theoretical (Petukhov) Data', 'Fitted Line');

figure(8);
hold on
errorbar(Re(31:39), f_experimental(31:39),Py4, 's','LineWidth',1.5); % error bar
plot(Re(31:39),exp(y4),'Color', [0 0 1], 'lineWidth', 2); % Fitted Line
plot(Re(31:39),exp(y4)+Py_hat4,'--','Color', [0 0 0], 'lineWidth', 2);
plot(Re(31:39),exp(y4)-Py_hat4,'--','Color', [0 0 0], 'lineWidth', 2);

hold off
xlabel('log(Re)');
ylabel('log(fr) [Friction Factor]');
title('Moody Plot for Square Rough Pipe');
legend('Error Bar (Precision uncertainty of a measurement)', 'Fitted Line','Error Lines (Precision uncertainity of the curve fit)');



%%

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


