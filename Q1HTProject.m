clc
clear
%% Data
s = @(t) 59116.31 + 58.03 * t - 0.138 * t^2 + 1.102 * 10^-4 * t^3 - 3.751 * 10^-8 * t^4 + 4.683 * 10^-12 * t^5;
Timbient = 25;
h = 10;
R = 9 * 10^-3;
Ru = 2939;
k = 1.6;
cp = 2400;
a = 2.268 * 10^-7;
r1 = 4 * 10^-3;
r2 = 8 * 10^-3;
t1 = 900;
t2 = 1800;
t3 = 3600;
dr = 0.5 * 10^-3;
dt = [0.5 0.4567 0.4566 0.4565 0.4];
%% Plotting
figure('Name', 'Question 1', 'NumberTitle', 'off')
for z = 1:length(dt)
        tmax = 3600;
        m = round(R/dr);
        n = round(tmax/dt(z));
        A = zeros(n + 1, m + 1);
        K = (a * dt(z)) / dr^2;
        t = linspace(0, tmax, n + 1);
        redius = 0:dr:R;
        %% Explicit
        A(1,: ) = Timbient;
        for p = 2 : n + 1
            for i = 1 : m + 1
                if i == 1
                    A(p, i) = ((a * dt(z)) / k) * s((p - 1 - 1) * dt(z)) + (1 - 4 * K) * A(p - 1, i) + (4 * K) * A(p - 1, i + 1); 
                elseif 1 < i && i < m + 1
                    A(p, i) = ((a * dt(z)) / k) * s((p - 1 - 1) * dt(z)) + (K - K/(2 * (i - 1))) * A(p - 1, i - 1) + (1 - 2 * K) * A(p - 1, i) + (K + K/(2 * (i - 1))) * A(p - 1, i + 1);
                elseif i == m  + 1
                    A(p, i) = ((a * dt(z)) / k) * s((p - 1 - 1) * dt(z)) + (K + K/(2*(i-1))) * ((2 * h * dr * Timbient)/k) + (2 * K) * A(p - 1, i -1) + (1 - 2 * K + (K + K/(2*(i-1))) * ((-2 * h * dr)/k)) * A(p - 1, i);
                end
            end
        end
        subplot(5, 1, z)
        plot(t, A( : , 19))
        title(['dt = ', num2str(dt(z)), ' and K = ', num2str(K)])
        xlabel('Time(s)')
        ylabel('Temperature(°C) in last node')
end
%% Finding The critical value of K and dt
disp('The critical value of K and dt with error and try are 0.41414 and 0.4565s')

 