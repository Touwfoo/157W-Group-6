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


%Values from lab manual
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
R=18.82*10^-2; %Units: m
m=[1,2,3,4,5]*0.45359237; %Units: kg
M=8*0.45359237; %Units: kg
g=9.81; %Units: m/s^2
E=68.3*10^9; %Units: Pa
v=0.33; %Unitless
G=E/2/(1+v); %Units: Pa

%Calculations for Pure Bending Pure Torsion, and Axial Distribution
%Note: Negative signs were added to fix sign conventions. Why? IDK

Mb=m*g*(L/2-11/16*(L-xload)); %Units: N*m
epsilon_b=32*Mb*D/pi/E/(D^4-d^4); %Unitless
Mt=-m*g*R; %Units: N*m
epsilon_t=16*Mt*D/pi/G/(D^4-d^4); %Unitless

Mb_axial(4:5)=5*M*g*(L-x(4:5))/16; %Units: N*m
Mb_axial(1:3)=M*g*(L/2-11/16*(L-x(1:3))); %Units: N*m

Mb_axial2(4:5)=5*M*g*(L-x(4:5))/16; %Units: N*m
Mb_axial2(1:3)=M*g*(L/2-11/16*(L-x(1:3))); %Units: N*m

M_experimental=-avgStrainAxialDistribution*pi*E*(D^4-d^4)/32/D;

x_values=linspace(0,L,100);
M_theoretical=zeros(size(x_values));
for i = 1:100
     if L-x_values(i) < L/2
        M_theoretical(i)=5*M*g*(L-x_values(i))/16; %Units: N*m
     else
        M_theoretical(i)=M*g*(L/2-11/16*(L-x_values(i))); %Units: N*m
     end
end
Epsilon_theoretical = -32*M_theoretical*D/pi/E/(D^4-d^4);


%Calculations for Cross-Talk Coefficients
k_values = zeros(3,4);
for n=1:3
    et1=avgStrainFixedBending(n,1);
    eb1=avgStrainFixedBending(n,2);
    et2=avgStrainFixedTorsion(n,1);
    eb2=avgStrainFixedTorsion(n,2);
    etprime1=avgStrainPureTorsion(n+1);
    ebprime1=avgStrainPureBending(1);
    etprime2=avgStrainPureTorsion(1);
    ebprime2=avgStrainPureBending(n+1);
    e_vector=[et1; eb1; et2; eb2];
    matrix=[etprime1, ebprime1, 0, 0;
            0, 0, etprime1, ebprime1;
            etprime2, ebprime2, 0, 0;
            0, 0, etprime2, ebprime2];
    k_vector=matrix\e_vector;
    k_values(n, :)=k_vector;
end
%Rows separated by weight, columns in order of ktt, ktb, kbt, kbb
k_values

%Error Analysis
%Given constants
deltam = 0.005; %Units: kg
deltaL = 0.0002*10^-2; %Units: m
deltaD = 0.0002*10^-2; %Units: m
deltad = 0.0002*10^-2; %Units: m
% lines of best fit
pureBendingBestFit = polyfit(Mb, avgStrainPureBending, 1);
pureTorsionBestFit = polyfit(Mt, avgStrainPureTorsion, 1);
MexpAxialBestFit1 = polyfit(x(1:3), M_experimental(1:3), 1);
MexpAxialBestFit2 = polyfit(x(3:5), M_experimental(3:5), 1);
StrainExpAxialBestFit1 = polyfit(x(1:3), avgStrainAxialDistribution(1:3), 1);
StrainExpAxialBestFit2 = polyfit(x(3:5), avgStrainAxialDistribution(3:5), 1);


%Plot for Pure Bending
figure(1);
hold on
plot(Mb, avgStrainPureBending, 'x','Color', [0 1 0], 'lineWidth', 2); %Experimental
plot(Mb, epsilon_b, 'o','Color', [1 0 0], 'lineWidth', 2); %Theoretical
plot(Mb, Mb*pureBendingBestFit(1) + pureBendingBestFit(2), '--', 'Color', [0 0 1], 'lineWidth', 2) %Line of best fit

hold off
grid on
xlabel('Moment [N*m]');
ylabel('Strain [-]');
title('Pure Bending');
legend('Experimental Data', 'Theoretical Data', 'Line of Best Fit');

%Plot for Pure Torsion
figure(2);
hold on
plot(Mt, avgStrainPureTorsion, 'x','Color', [0 1 0], 'lineWidth', 2); %Experimental
plot(Mt, epsilon_t, 'o','Color', [1 0 0], 'lineWidth', 2); %Theoretical
plot(Mt, Mt*pureTorsionBestFit(1) + pureTorsionBestFit(2), '--', 'Color', [0 0 1], 'lineWidth', 2); %Line of best fit


hold off
grid on
xlabel('Torque [N*m]');
ylabel('Strain [-]');
title('Pure Torsion');
legend('Experimental Data', 'Theoretical Data', 'Line of Best Fit');

%Plot for Axial Distribution
figure(3);
hold on
plot(x, M_experimental, 'x','Color', [0 1 0], 'lineWidth', 2); %Experimental
plot(x_values, M_theoretical, '-','Color', [1 0 0], 'lineWidth', 2); %Theoretical
plot(x_values(1:round(length(x_values)/2)), x_values(1:round(length(x_values)/2))*MexpAxialBestFit1(1) + MexpAxialBestFit1(2), '--', 'Color', [0 0 1], 'lineWidth', 2) %Line of Best Fit 1
plot(x_values(round(length(x_values)/3):end), x_values(round(length(x_values)/3):end)*MexpAxialBestFit2(1) + MexpAxialBestFit2(2), '--', 'Color', [0 0 1], 'lineWidth', 2) %Line of Best Fit 2
yline(0);

hold off
grid on
xlabel('Position [m]');
ylabel('Moment [N*m]');
title('Axial Distribution Moment Graph');
legend('Experimental Data', 'Theoretical Data', 'Line of Best Fit');

figure(4);
hold on
plot(x, avgStrainAxialDistribution, 'x','Color', [0 1 0], 'lineWidth', 2); %Experimental
plot(x_values, Epsilon_theoretical, '-','Color', [1 0 0], 'lineWidth', 2); %Theoretical
plot(x_values(1:round(length(x_values)/2)), x_values(1:round(length(x_values)/2))*StrainExpAxialBestFit1(1) + StrainExpAxialBestFit1(2), '--', 'Color', [0 0 1], 'lineWidth', 2) %Line of Best Fit 1
plot(x_values(round(length(x_values)/3):end), x_values(round(length(x_values)/3):end)*StrainExpAxialBestFit2(1) + StrainExpAxialBestFit2(2), '--', 'Color', [0 0 1], 'lineWidth', 2) %Line of Best Fit 2

yline(0);

hold off
grid on
xlabel('Position [m]');
ylabel('Strain [-]');
title('Axial Distribution Strain Graph');
legend('Experimental Data', 'Theoretical Data', 'Line of Best Fit');