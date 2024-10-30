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
hold on
for a = 1:7
    plot(time_nplc(:,a), temp_nplc(:,a));
end
hold off
legend('1', '2', '3', '4', '5', '7', '10')

figure(2)
hold on
for a = 1:3
    plot(time_normal(:,a), temp_normal(:,a));
end
hold off

figure(3)
hold on
for a = 1:3
    plot(time_clean(:,a), temp_clean(:,a));
end
hold off

figure(4)
hold on
for a = 1:3
    plot(time_vas(:,a), temp_vas(:,a));
end
hold off

%% plot derivatives
figure(5)
hold on
for a = 1:7
    plot(time_nplc(:,a), Dtdt_nplc(:,a));
end
hold off
legend('1', '2', '3', '4', '5', '7', '10')

figure(6)
hold on
for a = 1:3
    plot(time_normal(:,a), Dtdt_normal(:,a));
end
hold off

figure(7)
hold on
for a = 1:3
    plot(time_clean(:,a), Dtdt_clean(:,a));
end
hold off

figure(8)
hold on
for a = 1:3
    plot(time_vas(:,a), Dtdt_vas(:,a));
end
hold off

%% plot temp vs. heat transfer rates

figure(9)
for a = 1:7
    loglog(temp_nplc(:,a),-Dtdt_nplc(:,a));
    hold on
end
hold off

figure(10)
for a = 1:3
    loglog(temp_normal(:,a),-Dtdt_normal(:,a));
    hold on
end
hold off

figure(11)
for a = 1:3
    loglog(temp_clean(:,a),-Dtdt_clean(:,a));
    hold on
end
hold off

figure(12)
for a = 1:3
    loglog(temp_vas(:,a),-Dtdt_vas(:,a));
    hold on
end
hold off

%% perform moving average on the data, plot the data

Dtdt_moveAvg_nplc = getMoveAvg(Dtdt_nplc, 16);
Dtdt_moveAvg_normal = getMoveAvg(Dtdt_normal, 16);
Dtdt_moveAvg_clean = getMoveAvg(Dtdt_clean, 16);
Dtdt_moveAvg_vas = getMoveAvg(Dtdt_vas, 16);

figure(13)
for a = 1:7
    loglog(temp_nplc(:,a),-Dtdt_moveAvg_nplc(:,a));
    hold on
end
hold off
ylim([1e-1, 1e2]);

figure(14)
for a = 1:3
    loglog(temp_normal(:,a),-Dtdt_moveAvg_normal(:,a));
    hold on
end
hold off
ylim([1e-1, 1e2]);

figure(15)
for a = 1:3
    loglog(temp_clean(:,a),-Dtdt_moveAvg_clean(:,a));
    hold on
end
hold off
ylim([1e-1, 1e2]);

figure(16)
for a = 1:3
    loglog(temp_vas(:,a),-Dtdt_moveAvg_vas(:,a));
    hold on
end
hold off
ylim([1e-1, 1e2]);

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