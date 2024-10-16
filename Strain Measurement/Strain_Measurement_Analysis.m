%{
MECH&AE 157W - Basic Mechanical and Aerospace Engineering Laboratory with Writing
Instructor: Prof. Yongie Hu
Teaching Assistants: Zihao Qin, Qiyu Xing, Benjamin Heronimus
Fall 2024 – UCLA
Students: Alex Lie, Andrew Tan, Anli Liu, Umer Badae, Ian Lee
Lab Group: 6
Lab: Strain Measurement
%}

%Beginning of code
clear
close
clc
format long

%{
Row Information
Pure Bending: 1-25
Pure Torsion: 26-50
Fixed Bending and Variable Torsion: 51-65
Variable Bending and Fixed Torsion: 66-80
Axial Distribution: 81-85
%}

%Importing data and converting to appropriate units
data=readtable("Strain_Measurement_Data.xlsx");
BendingWeight=data.BendingWeight*0.45359237; %Units: kg
TorsionWeight=data.TorsionWeight*0.45359237; %Units: kg
TorsionReading=data.TorsionReading*10^-6; %Unitless
BendingReading1=data.BendingReading1*10^-6; %Unitless
BendingReading2=data.BendingReading2*10^-6;	%Unitless
BendingReading3=data.BendingReading3*10^-6;	%Unitless
BendingReading4=data.BendingReading4*10^-6;	%Unitless
BendingReading5=data.BendingReading5*10^-6; %Unitless

%In order of 1lb, 2lb, 3lb, 4lb, 5lb
avgStrainPureBending = [mean(BendingReading1(1:5)); 
                        mean(BendingReading1(6:10)); 
                        mean(BendingReading1(11:15)); 
                        mean(BendingReading1(16:20)); 
                        mean(BendingReading1(21:25))];

stdStrainPureBending = [std(BendingReading1(1:5)); 
                        std(BendingReading1(6:10)); 
                        std(BendingReading1(11:15)); 
                        std(BendingReading1(16:20)); 
                        std(BendingReading1(21:25))];

avgStrainPureTorsion = [mean(TorsionReading(26:30)); 
                        mean(TorsionReading(31:35)); 
                        mean(TorsionReading(36:40)); 
                        mean(TorsionReading(41:45)); 
                        mean(TorsionReading(46:50))];

stdStrainPureTorsion = [std(TorsionReading(26:30)); 
                        std(TorsionReading(31:35)); 
                        std(TorsionReading(36:40)); 
                        std(TorsionReading(41:45)); 
                        std(TorsionReading(46:50))];

%Columns in order of first torsion then bending, rows in order of 2lb, 3lb, 4lb
avgStrainFixedBending = [mean(TorsionReading(51:55)), mean(BendingReading1(51:55));
                        mean(TorsionReading(56:60)), mean(BendingReading1(56:60)); 
                        mean(TorsionReading(61:65)), mean(BendingReading1(61:65))];

stdStrainFixedBending = [std(TorsionReading(51:55)), std(BendingReading1(51:55));
                        std(TorsionReading(56:60)), std(BendingReading1(56:60)); 
                        std(TorsionReading(61:65)), std(BendingReading1(61:65))];

avgStrainFixedTorsion = [mean(TorsionReading(66:70)), mean(BendingReading1(66:70));
                        mean(TorsionReading(71:75)), mean(BendingReading1(71:75)); 
                        mean(TorsionReading(76:80)), mean(BendingReading1(76:80))];

stdStrainFixedTorsion = [std(TorsionReading(66:70)), std(BendingReading1(66:70));
                        std(TorsionReading(71:75)), std(BendingReading1(71:75)); 
                        std(TorsionReading(76:80)), std(BendingReading1(76:80))];

%{
%In order of x1, x2, x3, x4, x5
avgStrainAxialDistribution = [mean(BendingReading1(81:85)); 
                            mean(BendingReading2(81:85)); 
                            mean(BendingReading3(81:85)); 
                            mean(BendingReading4(81:85)); 
                            mean(BendingReading5(81:85))];

stdStrainAxialDistribution = [std(BendingReading1(81:85)); 
                            std(BendingReading2(81:85)); 
                            std(BendingReading3(81:85)); 
                            std(BendingReading4(81:85)); 
                            std(BendingReading5(81:85))];
%}

%Values from lab manual (not sure if accurate)
x1=1.5*0.0254; %Units: m
x2=3.5*0.0254; %Units: m
x3=12.375*0.0254; %Units: m
x4=24.375*0.0254; %Units: m
x5=30.35*0.0254; %Units: m
xload=15.313*0.0254; %Units: m
x=[x1; x2; x3; x4; x5];

%Aluminum Pipe Parameters
L=80.2*10^-2; %Units: m
D=1.905*10^-2; %Units: m
d=1.65*10^-2; %Units: m
R=18.882*10^-2; %Units: m   <---------------  IT MIGHT BE 18.82 cm???
Weight=[1,2,3,4,5]*0.45359237; %Units: kg
F=9.81*Weight; %Units: N
E=68.3*10^9; %Units: Pa
v=0.33; %Unitless
G=E/2/(1+v); %Units: Pa
I=(D^4-d^4)*pi/64; %Units: m^4
J=(D^4-d^4)*pi/32; %Units: m^4

%Calculations
Torque=F*R; %Units: N*m
Moment=F*L; %Units: N*m

%Plot for Pure Bending
figure(1);
hold on
plot(Moment, avgStrainPureBending, 'x','Color', [0 1 0], 'lineWidth', 2); %Experimental
plot(Moment, avgStrainPureBending, 'o','Color', [1 0 0], 'lineWidth', 2); %Theoretical

hold off
xlabel('Moment [N*m]');
ylabel('Strain [—]');
title('Pure Bending');
legend('Experimental Data', 'Theoretical Data');

%Plot for Pure Torsion
figure(2);
hold on
plot(Torque, avgStrainPureTorsion, 'x','Color', [0 1 0], 'lineWidth', 2); %Experimental
plot(Torque, avgStrainPureTorsion, 'o','Color', [1 0 0], 'lineWidth', 2); %Theoretical

hold off
xlabel('Torque [N*m]');
ylabel('Strain [—]');
title('Pure Torsion');
legend('Experimental Data', 'Theoretical Data');

%{
%Plot for Axial Distribution
figure(3);
hold on
plot(x, avgStrainAxialDistribution, 'x','Color', [0 1 0], 'lineWidth', 2); %Experimental
plot(x, avgStrainAxialDistribution, 'o','Color', [1 0 0], 'lineWidth', 2); %Theoretical

hold off
xlabel('Position [m]');
ylabel('Strain [—]');
title('Pure Torsion');
legend('Experimental Data', 'Theoretical Data');
%}

%Will work on theoretical stuff later