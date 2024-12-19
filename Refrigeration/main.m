%{
MECH&AE 157W - Basic Mechanical and Aerospace Engineering Laboratory with Writing
Instructor: Prof. Yongie Hu
Teaching Assistants: Zihao Qin, Qiyu Xing, Benjamin Heronimus
Fall 2024 – UCLA
Students: Alex Lie, Andrew Tan, Anli Liu, Umer Badae, Ian Lee
Lab Group: 6
Lab: PRefrigeration
%}

%Beginning of code
clear all
close all
clc

%Importing data
data = readmatrix("pure_data.csv");
%{
columns: P1-P4, T1-T8, T10, mdot, I1-I3, V
rows:
    pressure controlled: 2psig
    pressure controlled: 7psig
    pressure controlled: 15psig
    pressure controlled: 30psig
    thermal controlled: high evaporator, high condensor
    thermal controlled: high, low
    capillary tube: high, high
    capillary tube: high, low
    capillary tube: low, high
%}

