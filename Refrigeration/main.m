%{
MECH&AE 157W - Basic Mechanical and Aerospace Engineering Laboratory with Writing
Instructor: Prof. Yongie Hu
Teaching Assistants: Zihao Qin, Qiyu Xing, Benjamin Heronimus
Fall 2024 – UCLA
Students: Alex Lie, Andrew Tan, Anli Liu, Umer Badae, Ian Lee
Lab Group: 6
Lab: PRefrigeration
%}

%% set up
clear all
close all
clc

%% Import data
data = readmatrix("pure_data.csv");
% change gauge pressure to absolute pressure
data(:,2:3) = data(:,2:3) + 14.7;
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

R12properties = readmatrix("R12_properties.csv");
R12properties(:,2) = R12properties(:,2) * 14.5038; % unit conversion from bara to psia
% columns: T[C], P_s[psia], v_i, h_f, h_g, S_f, S_g, h@15Cover, S@15Cover, h@30Cover, S@30Cover

%% main code
for n = 1:1:9
    Trial1 = processTrial(n, data, R12properties);
    plotTrial(n, Trial1, R12properties)
end

%% functions
function out = approxSatProp(T, properties)
    % input temperature in celcius
    % output T[C] P_s[psia], v_i, h_f, h_g, S_f, S_g
    if T < -100
        out = properties(1,:);
        % approximate at -100 *C
    elseif T > 112
        out = properties(end,:);
        % approximate at 112 *C
    else
        % interpolation
        n = 1;
        while T > properties(n,1)
            n = n + 1;
        end
        X_0 = properties(n-1,:); % properties at temperature below
        X_1 = properties(n, :); % properties at temperature above
        out = X_0 + (X_1 - X_0).*(T-X_0(1))./(X_1(1)-X_0(1)); % linear approximation
    end
end

% approximating the properties of superheated vapor
function out = approxPropByP(P, properties)
    % input pressure in psia
    % output T[C] P_s[psia], v_i, h_f, h_g, S_f, S_g
     temp_prop = properties(1:0).*0;
    % basically the same approximation for sat except with pressure
    if P < properties(1, 2)
        temp_prop = properties(1,:);
        % approximate at -100 *C
    elseif P > properties(end, 2)
        temp_prop = properties(end,:);
        % approximate at 112 *C
    else
        % interpolation
        n = 1;
        while P > properties(n,2)
            n = n + 1;
        end
        X_0 = properties(n-1,:); % properties at temperature below
        X_1 = properties(n, :); % properties at temperature above
        temp_prop = X_0 + (X_1 - X_0).*(P-X_0(2))./(X_1(2)-X_0(2)); % linear approximation
    end
    out = temp_prop;
end

function out = approxSHPropByT(T, P, properties)
    % input Temperature in Celcius, Pressure is psia
    % output T, P, h, S
    temp_prop = approxPropByP(P, properties);
    % interpolate again based on temperature
    if T < temp_prop(1)
        % that's not a superheated vapor, approximate as sat vapor
        out = [T, P, temp_prop(5), temp_prop(7)];
        return
    elseif T - temp_prop(1) <= 15
        % between 0 and 15 *C over the sat temp
        hS_0 = [temp_prop(5), temp_prop(7)];
        hS_1 = [temp_prop(8), temp_prop(9)];
        temp_hS = hS_0 + (hS_1 - hS_0).*(T-temp_prop(1))./15;
    else
        % between 15 and 30 *C over the sat temp
        % if more than 30*C over, assume linear trend
        hS_0 = [temp_prop(8), temp_prop(9)];
        hS_1 = [temp_prop(10), temp_prop(11)];
        temp_hS = hS_0 + (hS_1 - hS_0).*(T-(temp_prop(1) + 15))./15;
    end
    h = temp_hS(1);
    S = temp_hS(2);
    out = [T, P, h, S];
end

function out = approxSHPropByS(S, P, properties)
    % input Temperature in Celcius, Pressure is psia
    % output T, P, h, S
    temp_prop = approxPropByP(P, properties);
    % interpolate again based on entropy
    if S < temp_prop(7)
        % that's not a superheated vapor, approximate as sat vapor
        out = [temp_prop(1), P, temp_prop(5), S];
        return
    elseif S <= temp_prop(9)
        % between entoropy for 0 and 15 *C over the sat temp
        Th_0 = [temp_prop(1), temp_prop(5)];
        Th_1 = [temp_prop(1)+15, temp_prop(8)];
        temp_Th = Th_0 + (Th_1 - Th_0).*(S-temp_prop(7))./(temp_prop(9)-temp_prop(7));
    else
        % between entoropy 15 and 30 *C over the sat temp
        % if entorpy would be more than 30*C over, assume linear trend
        Th_0 = [temp_prop(1)+15, temp_prop(8)];
        Th_1 = [temp_prop(1)+30, temp_prop(10)];
        temp_Th = Th_0 + (Th_1 - Th_0).*(S-temp_prop(9))./(temp_prop(11)-temp_prop(9));
    end
    T = temp_Th(1);
    h = temp_Th(2);
    out = [T, P, h, S];
end

function out = approxSHPropByh(h, P, properties)
    % input Temperature in Celcius, Pressure is psia
    % output T, P, h, S
    temp_prop = approxPropByP(P, properties);
    % interpolate again based on entropy
    if h < temp_prop(5)
        % that's not a superheated vapor, approximate as sat vapor
        out = [temp_prop(1), P, h, temp_prop(7)];
        return
    elseif h <= temp_prop(8)
        % between enthalpy for 0 and 15 *C over the sat temp
        TS_0 = [temp_prop(1), temp_prop(7)];
        TS_1 = [temp_prop(1)+15, temp_prop(9)];
        temp_TS = TS_0 + (TS_1 - TS_0).*(h-temp_prop(5))./(temp_prop(8)-temp_prop(5));
    else
        % between enthalpy 15 and 30 *C over the sat temp
        % if enthalpy would be more than 30*C over, assume linear trend
        TS_0 = [temp_prop(1)+15, temp_prop(9)];
        TS_1 = [temp_prop(1)+30, temp_prop(11)];
        temp_TS = TS_0 + (TS_1 - TS_0).*(h-temp_prop(8))./(temp_prop(10)-temp_prop(8));
    end
    T = temp_TS(1);
    S = temp_TS(2);
    out = [T, P, h, S];
end

function out = processTrial(Trial, pure_data, properties)
    % extract only the relevant data
    data = pure_data(Trial,:);
    % Reminder: columns: P1-P4, T1-T8, T10, mdot, I1-I3, V

    % Point 1: after condenser, assume pure saturated liquid
    T1 = data(5);
    Pt1_properties = approxSatProp(T1, properties);
    % extract the enthalpy and entropy
    P1 = Pt1_properties(2);
    h1 = Pt1_properties(4);
    S1 = Pt1_properties(6);

    % Point 2: after the expander, assume isenthalpic
    T2 = data(6);
    Pt2_properties = approxSatProp(T2, properties);
    P2 = Pt2_properties(2);
    h2 = h1; % isenthalpic
    hf2 = Pt2_properties(4);
    hg2 = Pt2_properties(5);
    % caclulate quality
    x = (h2 - hf2)/(hg2 - hf2);
    Sf2 = Pt2_properties(6);
    Sg2 = Pt2_properties(7);
    S2 = x*Sg2 + (1 - x)*Sf2;

    % Point 3: after the evaporator
    T3 = data(7);
    P3 = data(3);
    % interpolate the superheated properties
    Pt3_properties = approxSHPropByT(T3, P3, properties);
    h3 = Pt3_properties(3);
    S3 = Pt3_properties(4);

    % Point 4: (experimental) after the compressor
    T4_exp = data(8);
    P4_exp = data(4);
    % interpolate the superheated properties
    Pt4_properties_exp = approxSHPropByT(T4_exp, P4_exp, properties);
    h4_exp = Pt4_properties_exp(3);
    S4_exp = Pt4_properties_exp(4);

    % Point 4: (ideal)
    S4_ideal = S3;
    P4_ideal = data(4);
    % interpolation based on S
    Pt4_properties_ideal = approxSHPropByS(S4_ideal, P4_ideal, properties);
    T4_ideal = Pt4_properties_ideal(1);
    h4_ideal = Pt4_properties_ideal(3);

    % Point 4: (real, no heat loss)
    V = data(18);
    I = data(16);
    m_dot = data(14) * 60 /0.453592;
    % unit check
    % V*A/ (lb/min) vs KJ/kg
    % * (1lb/0.453592kg) * (60s/min)
    h4_real = V*I/m_dot + h3;
    P4_real = data(4);
    % interpolation based on h
    Pt4_properties_real = approxSHPropByh(h4_real, P4_real, properties);
    T4_real = Pt4_properties_real(1);
    S4_real = Pt4_properties_real(4);

    out = [ T1, P1, h1, S1;
            T2, P2, h2, S2;
            T3, P3, h3, S3;
            T4_exp, P4_exp, h4_exp, S4_exp;
            T4_ideal, P4_ideal, h4_ideal, S4_ideal;
            T4_real, P4_real, h4_real, S4_real;];
end

function plotTrial(n, points, properties)
    figure(n) %T-S diagram
    plot(points(:,4), points(:,1), 'x', 'color', 'r', 'linewidth', 2)
    hold on
    plot(properties(:,6), properties(:,1), 'k', 'linewidth', 2)
    plot(properties(:,7), properties(:,1), 'k', 'linewidth', 2)
    hold off
end