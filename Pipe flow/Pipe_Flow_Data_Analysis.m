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
IDK1 = data.Pressure1;
PipeNumber = data.PipeNumber;
Temperature = data.Temperature; %Units: Kelvin
Frequency = data.Frequency; %Units: Hz????
IDK2 = data.Pressure2;

%Pipe Inner Diamaters (in order of 1, 2, 3, 4)
PipeDiameter = [8.15; 11.8; 17.5; 9.93]; %Units: mm
PipeCrossSectionalArea= pi*PipeDiameter.^2/4; %Units: mm^2

%Pipe Lengths (in order of 1, 2, 3, 4)
PipeLength = [0.781; 1.143; 1.98; 0.337]; %Units: m

