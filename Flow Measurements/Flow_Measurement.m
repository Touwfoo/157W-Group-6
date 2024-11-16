clc
close
clear

T = readtable('Flow_Measurement_Data.csv');

rotameter_percentage = T.Rotameter_OfFlow;
ambient_pressure = T.AmbientPressure_Pa_;
sonic_gauge_pressure = T.SonicGaugePressue_Pa_;
temp = T.Temperature__C_;
closed_position = T.ClosedPosition_in_;
current_position = T.CurrentPosition_in_;

oriface_plate_absolute_pressure = T.OrifacePlate_AbsolutePressure_Pa_;
orifice_plate_DP = [T.OrificePlateDP0_Pa_, T.OrificePlateDP1_Pa_, T.OrificePlateDP2_Pa_, T.OrificePlateDP3_Pa_, T.OrificePlateDP4_Pa_, T.OrificePlateDP5_Pa_];

turbine_absolute_pressure = T.Turbine_AbsolutePressure_Pa_;
turbine_freq = T.Count_Hz_;

laminar_absolute_pressure = T.Laminar_AbsolutePressure_Pa_;
laminar_pressure_drop = T.PressureDrop_in_;

venturi_absolute_pressure = T.Venturi_AbsolutePressure_Pa_ * 1.01362;
venturi_DP = [T.VenturiDP0_Pa_, T.VenturiDP1_Pa_, T.VenturiDP2_Pa_, T.VenturiDP3_Pa_, T.VenturiDP4_Pa_, T.VenturiDP5_Pa_];

%% venturi meter

mu = 186.1 / 1e7;
rho = 1.1614; % kg / m^3
venturi_D = 0.0520446; % m
venturi_d = [0.0258064; vent_dia(0.5 * venturi_D); vent_dia(venturi_D); vent_dia(1.5 * venturi_D); vent_dia(2 * venturi_D); 0.0520446]; % m
g = 9.81;
gamma = 1.4;
venturi_A = [pi * (venturi_D / 2).^2; pi * (venturi_d ./ 2).^2];
% ideal mass flow rate
m_dot_ideal = venturi_A(2) * sqrt((2 * g * rho * venturi_DP(:,2)) / (1 - (venturi_A(2) / venturi_A(1)))^2);
% reynolds number at inlet #
venturi_Re = (m_dot_ideal * venturi_D) / (venturi_A(1) * mu);
venturi_p1 = venturi_absolute_pressure;
venturi_p2 = venturi_absolute_pressure + venturi_DP(:,1);
% compressibility factor
venturi_Y = sqrt((venturi_p2 ./ venturi_p1).^(2 ./ gamma) .* (gamma ./ (gamma - 1)) .* ((1 - (venturi_p2 ./ venturi_p1).^((gamma - 1) ./ (gamma))) ./ (1 - (venturi_p2 ./ venturi_p1))) .* ((1 - (venturi_A(2) ./ venturi_A(1)).^2) ./ (1 - (venturi_A(2) ./ venturi_A(1)).^2 .* (venturi_p2 ./ venturi_p1).^(2 ./ gamma))));
% from discharge chart in manual, eyeballed from Re at throat:
C = [0;0;0;0;0;0.96;0.97;0.97;0.97;0.97;0.97];
% actual mass flow rate
m_dot_venturi = m_dot_ideal .* C .* venturi_Y;

% pressure recovery
venturi_pressure_recovered = (venturi_absolute_pressure + venturi_DP) ./ venturi_absolute_pressure;
tap_points = [0, 0.5 * venturi_D, 1 * venturi_D, 1.5 * venturi_D, 2 * venturi_D, 0.1450086];
scatter(tap_points, venturi_pressure_recovered(2,:));

%% orifice plate



%% functions

function res = vent_dia(d)
    venturi_angle = 7 * (pi / 180);
    res = 0.0258064 + 2 * tan(venturi_angle) * (d - 0.0130048);
end