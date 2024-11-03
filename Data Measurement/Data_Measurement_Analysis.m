clc
close
clear

T = readtable('DataMeasurementLab.csv');
% temp data
time_normal = [T.time_normal_1, T.time_normal_2, T.time_normal_3];
time_clean = [T.time_clean_1, T.time_clean_2, T.time_clean_3];
time_vas = [T.time_vas_1, T.time_vas_2, T.time_vas_3];
time_nplc = [T.time_normal_1, T.time_nplc2, T.time_nplc3, T.time_nplc4, T.time_nplc5, T.time_nplc7, T.time_nplc10];
% temp data
temp_normal = [T.temp_normal_1, T.temp_normal_2, T.temp_normal_3] + 273.15;
temp_clean = [T.temp_clean_1, T.temp_clean_2, T.temp_clean_3] + 273.15;
temp_vas = [T.temp_vas_1, T.temp_vas_2, T.temp_vas_3] + 273.15;
temp_nplc = [T.temp_normal_2, T.temp_nplc2, T.temp_nplc3, T.temp_nplc4, T.temp_nplc5, T.temp_nplc7, T.temp_nplc10] + 273.15;
% get derivatives
Dtdt_normal = getDerivative(time_normal,temp_normal);
Dtdt_clean = getDerivative(time_clean, temp_clean);
Dtdt_vas = getDerivative(time_vas, temp_vas);
Dtdt_nplc = getDerivative(time_nplc, temp_nplc);

%% plot temps
figure(1)
tiledlayout(2,2);
nexttile;
hold on
for a = 1:7
    plot(time_nplc(:,a), temp_nplc(:,a));
end
hold off
legend('1', '2', '3', '4', '5', '7', '10')
title('Time vs. Temperature for Varying nPLC');
xlabel('Time (s)');
ylabel('Temperature (C)');
grid on;

nexttile;
hold on
for a = 1:3
    plot(time_normal(:,a), temp_normal(:,a));
end
hold off
legend('Trial 1', 'Trial 2', 'Trial 3');
title('Time vs. Temperature for Normal Surface Conditions');
xlabel('Time (s)');
ylabel('Temperature (C)');
grid on;

nexttile;
hold on
for a = 1:3
    plot(time_clean(:,a), temp_clean(:,a));
end
hold off
legend('Trial 1', 'Trial 2', 'Trial 3');
title('Time vs. Temperature for Clean Surface Conditions');
xlabel('Time (s)');
ylabel('Temperature (C)');
grid on;

nexttile;
hold on
for a = 1:3
    plot(time_vas(:,a), temp_vas(:,a));
end
hold off
legend('Trial 1', 'Trial 2', 'Trial 3');
title('Time vs. Temperature for Vaseline Surface Conditions');
xlabel('Time (s)');
ylabel('Temperature (C)');
grid on;

%% plot derivatives
figure(2)
tiledlayout(2,2);

nexttile;
hold on
for a = 1:7
    plot(time_nplc(:,a), Dtdt_nplc(:,a));
end
hold off
legend('1', '2', '3', '4', '5', '7', '10')
title('Time vs. dTdt for Varying nPLC');
xlabel('Time (s)');
ylabel('Rate of Temperature Change (C/s)');
grid on;

nexttile;
hold on
for a = 1:3
    plot(time_normal(:,a), Dtdt_normal(:,a));
end
hold off
legend('Trial 1', 'Trial 2', 'Trial 3');
title('Time vs. dTdt for Normal Surface');
xlabel('Time (s)');
ylabel('Rate of Temperature Change (C/s)');
grid on;

nexttile;
hold on
for a = 1:3
    plot(time_clean(:,a), Dtdt_clean(:,a));
end
hold off
legend('Trial 1', 'Trial 2', 'Trial 3');
title('Time vs. dTdt for Clean Surface');
xlabel('Time (s)');
ylabel('Rate of Temperature Change (C/s)');
grid on;

nexttile;
hold on
for a = 1:3
    plot(time_vas(:,a), Dtdt_vas(:,a));
end
hold off
legend('Trial 1', 'Trial 2', 'Trial 3');
title('Time vs. dTdt for Vaseline Surface');
xlabel('Time (s)');
ylabel('Rate of Temperature Change (C/s)');
grid on;

%% 

D = 1.27 / 100; % diameter of the ball
rho = 8900; % density of copper
T_table = [70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300];
c_table = [171,203,229,251,271,287,301,313,323,332,339,346,351,356,365,372,377,382,386];
c_nplc = interp1(T_table, c_table,temp_nplc);
c_normal = interp1(T_table, c_table,temp_normal);
c_clean = interp1(T_table, c_table,temp_clean);
c_vas = interp1(T_table, c_table,temp_vas);
q_nplc = -(1/6) .* c_nplc .* rho .* D .* Dtdt_nplc;
q_normal = -(1/6) .* c_normal .* rho .* D .* Dtdt_normal;
q_clean = -(1/6) .* c_clean .* rho .* D .* Dtdt_clean;
q_vas = -(1/6) .* c_vas .* rho .* D .* Dtdt_vas;

T_sat = -195.8 + 273.15;
deltaT_nplc = temp_nplc - T_sat;
deltaT_normal = temp_normal - T_sat;
deltaT_clean = temp_clean - T_sat;
deltaT_vas = temp_vas - T_sat;

%% plot temp vs. heat transfer rates

figure(3)

tiledlayout(2,2);

nexttile;
for a = 1:7
    loglog(deltaT_nplc(:,a),q_nplc(:,a));
    hold on
end
hold off
ylim([2000, 200000]);
grid on;
title('Heat Flux vs. Temp, Varying nPLC');
xlabel('T-Tsat (K)');
ylabel('Heat Flux (W/m^2)');
legend('1', '2', '3', '4', '5', '7', '10');

nexttile;
for a = 1:3
    loglog(deltaT_normal(:,a),q_normal(:,a));
    hold on
end
hold off
ylim([2000, 200000]);
grid on;
title('Heat Flux vs. Temp, Normal Surface');
xlabel('T-Tsat (K)');
ylabel('Heat Flux (W/m^2)');
legend('Trial 1', 'Trial 2', 'Trial 3');

nexttile;
for a = 1:3
    loglog(deltaT_clean(:,a),q_clean(:,a));
    hold on
end
hold off
ylim([2000, 200000]);
grid on;
title('Heat Flux vs. Temp, Clean Surface');
xlabel('T-Tsat (K)');
ylabel('Heat Flux (W/m^2)');
legend('Trial 1', 'Trial 2', 'Trial 3');

nexttile;
for a = 1:3
    loglog(deltaT_vas(:,a),q_vas(:,a));
    hold on
end
hold off
ylim([2000, 200000]);
grid on;
title('Heat Flux vs. Temp, Vaseline Surface');
xlabel('T-Tsat (K)');
ylabel('Heat Flux (W/m^2)');
legend('Trial 1', 'Trial 2', 'Trial 3');

%% perform moving average on the data, plot the data

q_moveAvg_nplc = getMoveAvg(q_nplc, 16);
q_moveAvg_normal = getMoveAvg(q_normal, 16);
q_moveAvg_clean = getMoveAvg(q_clean, 16);
q_moveAvg_vas = getMoveAvg(q_vas, 16);

figure(4)
tiledlayout(2,2);

nexttile;
for a = 1:7
    loglog(deltaT_nplc(:,a),q_moveAvg_nplc(:,a));
    hold on
end
hold off
ylim([2000, 200000]);
grid on;
title('Heat Flux vs. Temp, Varying nPLC, w/ Moving Avg');
xlabel('T-Tsat (K)');
ylabel('Heat Flux (W/m^2)');
legend('1', '2', '3', '4', '5', '7', '10')

nexttile;
for a = 1:3
    loglog(deltaT_normal(:,a),q_moveAvg_normal(:,a));
    hold on
end
hold off
ylim([2000, 200000]);
grid on;
title('Heat Flux vs. Temp, Normal Surface, w/ Moving Avg');
xlabel('T-Tsat (K)');
ylabel('Heat Flux (W/m^2)');
legend('Trial 1', 'Trial 2', 'Trial 3');

nexttile;
for a = 1:3
    loglog(deltaT_clean(:,a),q_moveAvg_clean(:,a));
    hold on
end
hold off
ylim([2000, 200000]);
grid on;
title('Heat Flux vs. Temp, Clean Surface, w/ Moving Avg');
xlabel('T-Tsat (K)');
ylabel('Heat Flux (W/m^2)');
legend('Trial 1', 'Trial 2', 'Trial 3');

nexttile;
for a = 1:3
    loglog(deltaT_vas(:,a),q_moveAvg_vas(:,a));
    hold on
end
hold off
ylim([2000, 200000]);
grid on;
title('Heat Flux vs. Temp, Vaseline Surface, w/ Moving Avg');
xlabel('T-Tsat (K)');
ylabel('Heat Flux (W/m^2)');
legend('Trial 1', 'Trial 2', 'Trial 3');

%%
function moveAvg = getMoveAvg(data,k)
    Dtdt_moveAvg_ans = zeros(length(data(:,1)), length(data(1,:)));
    for a = 1:length(data(1,:))
        Dtdt_moveAvg_ans(:,a) = movmean(data(:,a), k);
    end
    moveAvg = Dtdt_moveAvg_ans;
end

function Dtdt = getDerivative(x,y)
    derivative = zeros(length(x(:,1)), length(x(1,:)));
    for b = 1:length(x(1,:))
        x_curr = x(:,b);
        y_curr = y(:,b);
        for a = 1:length(x_curr)
            if a < length(x_curr)
                derivative(a,b) = (y_curr(a+1) - y_curr(a)) / (x_curr(a + 1) - x_curr(a));
            else
                derivative(a,b) = derivative(a - 1,b);
            end
        end
    end
    Dtdt = derivative;
end
