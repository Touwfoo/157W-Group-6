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
% Pressue controlled
Trial1 = processTrial(1, data, R12properties); % 2psig
Trial2 = processTrial(2, data, R12properties); % 7psig
Trial3 = processTrial(3, data, R12properties); % 15psig
Trial4 = processTrial(4, data, R12properties); % 30psig
% Temperature controlled
Trial5 = processTrial(5, data, R12properties); % evap high, cond high
Trial6 = processTrial(6, data, R12properties); % high, low
% Capillary Tube
Trial7 = processTrial(7, data, R12properties); % high, high
Trial8 = processTrial(8, data, R12properties); % high, low
Trial9 = processTrial(9, data, R12properties); % low, high

% Plots
%{
hS and TS diagram for Pressure controlled at 15psig (2)
hS and TS diagram for Temperature controlled (4)
PV diagrams for Capillary Tube (3)
%}
%Pressure controlled
plot_hS(1, Trial3, R12properties)
xlabel('Entropy, S (kJ/kgK)')
ylabel('Enthalpy, h (kJ/kg)')
title('h-S Diagram for Pressure Controlled Expansion Valve','15psig')
legend('Data Points', 'Experimental', 'Ideal', 'Real', 'Saturation Dome')
xlim([0.2 1.2])
ylim([50 350])

plot_TS(2, Trial3, R12properties, data)
xlabel('Entropy, S (kJ/kgK)')
ylabel(['Temperature, T (' char(176) 'C)'])
title('T-S Diagram for Pressure Controlled Expansion Valve','15psig')
legend('Data Points', 'Experimental', 'Ideal', 'Real', 'Inlet Air Temperature', 'Saturation Dome')
xlim([0.2 1.2])
ylim([-20 250])

% Temperature controlled: high, high
plot_hS(3, Trial5, R12properties)
xlabel('Entropy, S (kJ/kgK)')
ylabel('Enthalpy, h (kJ/kg)')
title('h-S Diagram for Temperature Controlled Expansion Valve','Evaporator and Condenser Fans Set to High')
legend('Data Points', 'Experimental', 'Ideal', 'Real', 'Saturation Dome')
xlim([0.2 1])
ylim([50 350])

plot_TS(4, Trial5, R12properties, data)
xlabel('Entropy, S (kJ/kgK)')
ylabel(['Temperature, T (' char(176) 'C)'])
title('T-S Diagram for Temperature Controlled Expansion Valve','Evaporator and Condenser Fans Set to High')
legend('Data Points', 'Experimental', 'Ideal', 'Real', 'Inlet Air Temperature', 'Saturation Dome')
xlim([0.2 1])
ylim([-20 200])

% Temperature controlled: high, low
plot_hS(5, Trial6, R12properties)
xlabel('Entropy, S (kJ/kgK)')
ylabel('Enthalpy, h (kJ/kg)')
title('h-S Diagram for Temperature Controlled Expansion Valve','Evaporator Fan Set to High and Condenser Fan Set to Low')
legend('Data Points', 'Experimental', 'Ideal', 'Real','Saturation Dome')
xlim([0.2 1])
ylim([50 350])

plot_TS(6, Trial6, R12properties, data)
xlabel('Entropy, S (kJ/kgK)')
ylabel(['Temperature, T (' char(176) 'C)'])
title('T-S Diagram for Temperature Controlled Expansion Valve', 'Evaporator Fan Set to High and Condenser Fan Set to Low')
legend('Data Points', 'Experimental', 'Ideal', 'Real', 'Inlet Air Temperature',  'Saturation Dome')
xlim([0.2 1])
ylim([0 175])

% Capillary tube: all conditions
plot_PV(7, Trial7, R12properties)
xlabel('Specific Volume, v (m^3/kg)')
ylabel('Pressure, P (psia)')
title('P-v Diagram for Capillary Tube Expander','Evaporator and Condenser Fans Set to High')
legend('Data Points', 'Experimental', 'Ideal', 'Real', 'Saturation Dome')
xlim([-0.02 0.1])
ylim([30 200])

plot_PV(8, Trial8, R12properties)
xlabel('Specific Volume, v (m^3/kg)')
ylabel('Pressure, P (psia)')
title('P-v Diagram for Capillary Tube Expander','Evaporator Fan Set to High and Condenser Fan Set to Low')
legend('Data Points', 'Experimental', 'Ideal', 'Real', 'Saturation Dome')
xlim([-0.02 0.1])
ylim([30 250])

plot_PV(9, Trial9, R12properties)
xlabel('Specific Volume, v (m^3/kg)')
ylabel('Pressure, P (psia)')
title('P-v Diagram for Capillary Tube Expander','Evaporator Fan Set to Low and Condenser Fan Set to High')
legend('Data Points', 'Experimental', 'Ideal', 'Real', 'Saturation Dome')
xlim([-0.02 0.1])
ylim([30 200])

% calculate the COP, RC, RHR, Superheat of system
Calcs = zeros(9,7);
Calcs(1,:) = calculate_values(Trial1, data(1,:), R12properties);
Calcs(2,:) = calculate_values(Trial2, data(2,:), R12properties);
Calcs(3,:) = calculate_values(Trial3, data(3,:), R12properties);
Calcs(4,:) = calculate_values(Trial4, data(4,:), R12properties);
Calcs(5,:) = calculate_values(Trial5, data(5,:), R12properties);
Calcs(6,:) = calculate_values(Trial6, data(6,:), R12properties);
Calcs(7,:) = calculate_values(Trial7, data(7,:), R12properties);
Calcs(8,:) = calculate_values(Trial8, data(8,:), R12properties);
Calcs(9,:) = calculate_values(Trial9, data(9,:), R12properties);

% plot COP as func of inlet pressure
figure(10)
Inlet_Pressure = [2 7 15 30];
plot(Inlet_Pressure, Calcs(1:4,1), 'x', 'color', 'r', 'linewidth', 2)
hold on
plot(Inlet_Pressure, Calcs(1:4, 2), 'o', 'color', 'b', 'linewidth', 2)
COP_polyfit_ideal = polyfit(Inlet_Pressure, Calcs(1:4, 1), 1);
COP_polyfit_exp = polyfit(Inlet_Pressure, Calcs(1:4, 2), 1);
plot(Inlet_Pressure, COP_polyfit_ideal(1)*Inlet_Pressure + COP_polyfit_ideal(2), '--', 'color', 'r', 'linewidth', 1)
plot(Inlet_Pressure, COP_polyfit_exp(1)*Inlet_Pressure + COP_polyfit_exp(2), '--', 'color', 'b', 'linewidth', 1)
hold off
xlim([0 35])
xlabel('Inlet Pressure (psig)')
ylabel('COP')
title('COP of Pressure Regulated Expansion Valve at Varying Inlet Pressures')
legend('COP Ideal', 'COP Experimental', 'COP Ideal Trendline', 'COP Experimental Trendline')

% plot RC as a func of inlet pressure
figure(11)
plot(Inlet_Pressure, Calcs(1:4,3), 'x', 'color', 'r', 'linewidth', 2)
hold on
RC_polyfit = polyfit(Inlet_Pressure, Calcs(1:4, 3), 1);
plot(Inlet_Pressure, RC_polyfit(1)*Inlet_Pressure + RC_polyfit(2), '--', 'color', 'r', 'linewidth', 1)
hold off
xlim([0 35])
xlabel('Inlet Pressure (psig)')
ylabel('Refrigeration Capacity, RC')
title('Refrigeration Capacity of Pressure Regulated Expansion Valve at Varying Inlet Pressures')
legend('RC data', 'RC Trendline')

%% functions
function out = approxSatProp(T, properties)
    % input temperature in celcius
    % output T[C] P_s[psia], v_i, h_f, h_g, S_f, S_g
    if T < -100
        out = [properties(1,:), approxLiquidSV(-100)];
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
        X_0 = [properties(n-1,:), approxLiquidSV(properties(n-1,1))]; % properties at temperature below
        X_1 = [properties(n, :), approxLiquidSV(properties(n,1))]; % properties at temperature above
        out = X_0 + (X_1 - X_0).*(T-X_0(1))./(X_1(1)-X_0(1)); % linear approximation
    end
end

% approximating the properties of superheated vapor
function out = approxPropByP(P, properties)
    % input pressure in psia
    % output T[C] P_s[psia], v_i, h_f, h_g, S_f, S_g
     temp_prop = zeros(1, length(properties(1,:)) + 1);
    % basically the same approximation for sat except with pressure
    if P < properties(1, 2)
        temp_prop = [properties(1,:), approxLiquidSV(properties(1,1))];
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
        X_0 = [properties(n-1,:), approxLiquidSV(properties(n-1,1))]; % properties at temperature below
        X_1 = [properties(n, :), approxLiquidSV(properties(n,1))]; % properties at temperature above
        temp_prop = X_0 + (X_1 - X_0).*(P-X_0(2))./(X_1(2)-X_0(2)); % linear approximation
    end
    out = temp_prop;
end

function out = approxSHPropByT(T, P, properties)
    % input Temperature in Celcius, Pressure is psia
    % output T, P, h, S, vf, vg
    temp_prop = approxPropByP(P, properties);
    % interpolate again based on temperature
    if T < temp_prop(1)
        % that's not a superheated vapor, approximate as sat vapor
        out = [T, P, temp_prop(5), temp_prop(7), temp_prop(3)];
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
    v = (T+273)/(temp_prop(1)+273)*temp_prop(3)*0.995;
    out = [T, P, h, S, v];
end

function out = approxSHPropByS(S, P, properties)
    % input Temperature in Celcius, Pressure is psia
    % output T, P, h, S
    temp_prop = approxPropByP(P, properties);
    % interpolate again based on entropy
    if S < temp_prop(7)
        % that's not a superheated vapor, approximate as sat vapor
        out = [temp_prop(1), P, temp_prop(5), S, temp_prop(3)];
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
    v = (T+273)/(temp_prop(1)+273)*temp_prop(3)*0.995;
    out = [T, P, h, S, v];
end

function out = approxSHPropByh(h, P, properties)
    % input Temperature in Celcius, Pressure is psia
    % output T, P, h, S
    temp_prop = approxPropByP(P, properties);
    % interpolate again based on entropy
    if h < temp_prop(5)
        % that's not a superheated vapor, approximate as sat vapor
        out = [temp_prop(1), P, h, temp_prop(7), temp_prop(3)];
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
    v = (T+273)/(temp_prop(1)+273)*temp_prop(3)*0.995;
    out = [T, P, h, S, v];
end

% assume linear distribution for specific volume of liquid
function out = approxLiquidSV(T)
    % assume linear linear betweeen 250K and 300K
    tempT = T + 273;
    out = 1468^-1 + (1304^-1 - 1468^-1)*(tempT-250)/50;
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
    v1 = Pt1_properties(12);

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
    vf2 = Pt2_properties(12);
    vg2 = Pt2_properties(3);
    v2 = x*vg2 + (1 - x)*vf2;

    % Point 3: after the evaporator
    T3 = data(7);
    P3 = data(3);
    % interpolate the superheated properties
    Pt3_properties = approxSHPropByT(T3, P3, properties);
    h3 = Pt3_properties(3);
    S3 = Pt3_properties(4);
    v3 = Pt3_properties(5);

    % Point 4: (experimental) after the compressor
    T4_exp = data(9);
    P4_exp = data(4);
    % interpolate the superheated properties
    Pt4_properties_exp = approxSHPropByT(T4_exp, P4_exp, properties);
    h4_exp = Pt4_properties_exp(3);
    S4_exp = Pt4_properties_exp(4);
    v4_exp = Pt4_properties_exp(5);

    % Point 4: (ideal)
    S4_ideal = S3;
    P4_ideal = data(4);
    % interpolation based on S
    Pt4_properties_ideal = approxSHPropByS(S4_ideal, P4_ideal, properties);
    T4_ideal = Pt4_properties_ideal(1);
    h4_ideal = Pt4_properties_ideal(3);
    v4_ideal = Pt4_properties_ideal(5);

    % Point 4: (real, no heat loss)
    V = data(18);
    I = data(16);
    m_dot = data(14) /60 *0.453592;
    % unit check
    % V*A/ (lb/min) vs KJ/kg
    % * (1lb/0.453592kg) * (60s/min)
    h4_real = V*I/m_dot *1E-3 + h3;
    P4_real = data(4);
    % interpolation based on h
    Pt4_properties_real = approxSHPropByh(h4_real, P4_real, properties);
    T4_real = Pt4_properties_real(1);
    S4_real = Pt4_properties_real(4);
    v4_real = Pt4_properties_real(5);

    out = [ T1, P1, h1, S1, v1;
            T2, P2, h2, S2, v2;
            T3, P3, h3, S3, v3;
            T4_exp, P4_exp, h4_exp, S4_exp, v4_exp;
            T4_ideal, P4_ideal, h4_ideal, S4_ideal, v4_ideal;
            T4_real, P4_real, h4_real, S4_real, v4_real;];
end

function plot_quad(n_x, n_y, points)
    % create the dotted lines for each of the functions
    x = points(:,n_x);
    y = points(:,n_y);
    color = ['r', 'g', 'b'];
    for n = 1:1:3
        x_points = [x(1:3); x(3+n); x(1)];
        y_points = [y(1:3); y(3+n); y(1)];
        plot(x_points, y_points, '--', 'color', color(n), 'linewidth', 1)
    end
end

function plot_hS(n, points, properties)
    figure(n)
    plot(points(:,4), points(:,3), 'x', 'color', 'c', 'linewidth', 2)
    hold on
    plot_quad(4, 3, points)
    plot(properties(:,6), properties(:,4), 'k', 'linewidth', 1)
    plot(properties(:,7), properties(:,5), 'k', 'linewidth', 1)
    hold off
end

function plot_TS(n, points, properties, other_data)
    figure(n)
    plot(points(:,4), points(:,1), 'x', 'color', 'c', 'linewidth', 2)
    hold on
    plot_quad(4, 1, points)
    yline(other_data(12), '--', 'color', 'k', 'linewidth', 1)
    plot(properties(:,6), properties(:,1), 'k', 'linewidth', 1)
    plot(properties(:,7), properties(:,1), 'k', 'linewidth', 1)
    hold off
end

function plot_PV(n, points, properties)
    figure(n)
    plot(points(:,5), points(:,2), 'x', 'color', 'c', 'linewidth', 2)
    hold on
    plot_quad(5, 2, points)
    plot(approxLiquidSV(properties(:,1)), properties(:,2), 'k', 'linewidth', 1)
    plot(properties(:,3), properties(:,2), 'k', 'linewidth', 1)
    hold off
    xlim([0, max(points(:,5))*1.2])
end

function out = calculate_values(points, other_data, properties)
    m_dot = other_data(14) /60 *0.453592;
    %COP ideal, exp
    COP_ideal = (points(3,3) - points(2,3))/(points(5,3) - points(3,3));
    COP_exp = (points(3,3) - points(2,3))/(other_data(18)*other_data(16)/m_dot)*1E3;
    % REfrigeration Capacity
    RC = m_dot*(points(3,3) - points(2,3));
    % Rate of Heat Rejection
    RHR_exp = m_dot*(points(4,3) - points(1,1));
    RHR_ideal = m_dot*(points(5,3) - points(1,1));
    RHR_real = m_dot*(points(6,3) - points(1,1));
    % Superheating - amount of enthalpy over saturation - point 3
    temp = approxSatProp(points(3,1), properties); % to find sat vapor h_g
    Q_SH = m_dot*(points(3,3)- temp(4));
    out = [COP_ideal, COP_exp, RC, RHR_exp, RHR_ideal, RHR_real, Q_SH];
end