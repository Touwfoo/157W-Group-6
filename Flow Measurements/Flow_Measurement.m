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
sonic_h = 0.00254 .* (current_position - closed_position);

orifice_plate_absolute_pressure = T.OrifacePlate_AbsolutePressure_Pa_;
orifice_plate_DP = [T.OrificePlateDP0_Pa_, T.OrificePlateDP1_Pa_, T.OrificePlateDP2_Pa_, T.OrificePlateDP3_Pa_, T.OrificePlateDP5_Pa_];

turbine_absolute_pressure = T.Turbine_AbsolutePressure_Pa_;
turbine_freq = T.Count_Hz_;

laminar_absolute_pressure = T.Laminar_AbsolutePressure_Pa_;
laminar_pressure_drop = T.PressureDrop_in_;

venturi_absolute_pressure = T.Venturi_AbsolutePressure_Pa_ * 1.01362;
venturi_DP = [T.VenturiDP0_Pa_, T.VenturiDP1_Pa_, T.VenturiDP2_Pa_, T.VenturiDP3_Pa_, T.VenturiDP5_Pa_]; % skip 4

%% venturi meter

mu = 186.1 / 1e7;
rho = 1.1614; % kg / m^3
venturi_D = 0.0520446; % m
venturi_d = [0.0258064; vent_dia(0.5 * venturi_D); vent_dia(venturi_D); vent_dia(1.5 * venturi_D); 0.0520446]; % m
g = 9.81;
gamma = 1.4;
venturi_A = [pi * (venturi_D / 2).^2; pi * (venturi_d ./ 2).^2];
% ideal mass flow rate
m_dot_ideal = venturi_A(2) * sqrt((2 * g * rho * venturi_DP(:,2)) / (1 - (venturi_A(2) / venturi_A(1)))^2);
% reynolds number at inlet #
venturi_Re = (m_dot_ideal * venturi_D) / (venturi_A(1) * mu);
venturi_p1 = venturi_absolute_pressure;
venturi_p2 = venturi_absolute_pressure - venturi_DP(:,1);
% compressibility factor
venturi_Y = sqrt((venturi_p2 ./ venturi_p1).^(2 ./ gamma) .* (gamma ./ (gamma - 1)) .* ((1 - (venturi_p2 ./ venturi_p1).^((gamma - 1) ./ (gamma))) ./ (1 - (venturi_p2 ./ venturi_p1))) .* ((1 - (venturi_A(2) ./ venturi_A(1)).^2) ./ (1 - (venturi_A(2) ./ venturi_A(1)).^2 .* (venturi_p2 ./ venturi_p1).^(2 ./ gamma))));
% from discharge chart in manual, eyeballed from Re at throat:
C = [0;0;0;0;0;0.96;0.97;0.97;0.97;0.97;0.97];
% actual mass flow rate
m_dot_venturi = m_dot_ideal .* C .* venturi_Y;

% pressure recovery
venturi_pressure_recovered = (venturi_absolute_pressure - venturi_DP) ./ venturi_absolute_pressure;
tap_points = [0, 0.5 * venturi_D, 1 * venturi_D, 1.5 * venturi_D, 0.1450086];

figure(1);
hold on
for a = 1:11
    plot(tap_points, venturi_pressure_recovered(a,:), '-o');
end
hold off
title('Venturi Pressure Recovery');
xlabel('Tap position from throat (m)');
ylabel('Recovery %');
legend('0% flow', '10% flow', '20% flow', '30% flow', '40% flow', '50% flow', '60% flow', '70% flow', '80% flow', '90% flow', '96% flow')

%% orifice plate

orifice_d = 0.0254;
orifice_D = 0.0520446;
orifice_A = [2 * (orifice_D / 2)^2, 2 * (orifice_d / 2)^2];
orifice_m_dot_ideal = orifice_A(2) * sqrt((2 * rho * orifice_plate_DP(:,1)) / (1 - (orifice_A(2) / orifice_A(1))^2));
orifice_P = orifice_plate_absolute_pressure - orifice_plate_DP;
orifice_Y = 1 - (0.41 + 0.35 * (orifice_d / orifice_D)^4) * (1 - orifice_plate_absolute_pressure ./ orifice_P(:,1)) * (1 / gamma);
orifice_C = m_dot_venturi ./ (orifice_Y .* orifice_m_dot_ideal);

figure(2);
scatter(venturi_Re, orifice_C, 'filled');
title('Orifice Plate Dicharge Coefficient vs. Reynolds Number');
xlabel('Reynolds number');
ylabel('Discharge coefficient');

orifice_tap_points = orifice_D * [0.5, 1, 1.5, 2, 4];
orifice_pressure_recovered = (orifice_plate_absolute_pressure - orifice_plate_DP) ./ orifice_plate_absolute_pressure;

figure(3);
hold on
for a = 1:11
    plot(orifice_tap_points, orifice_pressure_recovered(a,:), '-o');
end
hold off
title('Orifice Plate Pressure Recovery vs. Location');
xlabel('Tap position from throat (m)');
ylabel('Recovery %');
legend('0% flow', '10% flow', '20% flow', '30% flow', '40% flow', '50% flow', '60% flow', '70% flow', '80% flow', '90% flow', '96% flow')

%% sonic nozzle

sonic_D = 0.0079248;
sonic_angle = 8.88 * (pi / 180);
sonic_A = (pi / 4) .* (sonic_D^2 - (sonic_D - 2 .* sonic_h .* tan(sonic_angle)).^2);
sonic_m_dot = sonic_A * (202438.7) * sqrt((2 * g) / (286.9 * 294.15)) * sqrt((gamma / (gamma + 1)) * (2 / (gamma + 1))^(2 / (gamma - 1)));
sonic_C = m_dot_venturi ./ sonic_m_dot;

figure(4);
scatter(venturi_Re, sonic_C, 'filled');
title('Sonic Nozzle Dicharge Coefficient vs. Reynolds Number');
xlabel('Reynolds number');
ylabel('Discharge coefficient');

%% laminar flow meter



%% functions

function res = vent_dia(d)
    venturi_angle = 7 * (pi / 180);
    res = 0.0258064 + 2 * tan(venturi_angle) * (d - 0.0130048);
end