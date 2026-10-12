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
dt = 0.25;
tmax = 3600;
m = round(R/dr);
n = round(tmax/dt);
A = zeros(n + 1, m + 1);
Z = zeros(n + 1, m + 1);
K = (a * dt) / dr^2;
t = 0:dt:tmax;
redius = 0:dr:R;
%% Explicit
A(1,: ) = Timbient;
for p = 2 : n + 1
    for i = 1 : m + 1
        if i == 1
            A(p, i) = ((a * dt) / k) * s((p - 1 - 1) * dt) + (1 - 4 * K) * A(p - 1, i) + (4 * K) * A(p - 1, i + 1);
        elseif 1 < i && i < m + 1
            A(p, i) = ((a * dt) / k) * s((p - 1 - 1) * dt) + (K - K/(2 * (i - 1))) * A(p - 1, i - 1) + (1 - 2 * K) * A(p - 1, i) + (K + K/(2 * (i - 1))) * A(p - 1, i + 1);
        elseif i == m  + 1
            A(p, i) = ((a * dt) / k) * s((p - 1 - 1) * dt) + (K + K/(2*(i-1))) * ((2 * h * dr * Timbient)/k) + (2 * K) * A(p - 1, i -1) + (1 - 2 * K + (K + K/(2*(i-1))) * ((-2 * h * dr)/k)) * A(p - 1, i);
        end
    end
end
%% Crank–Nicolson
Z(1,: ) = Timbient;
C = zeros(m + 1 , m + 1);
D = zeros(m + 1, 1);
for p = 2 : n + 1
    for i = 1 : m + 1
       if i == 1
            C(1 , 1) = 1 + 2 * K;
            C(1 , 2) = -2 * K;
            D(1, 1) = ((a * dt) / (2 * k)) * (s((p - 1 - 1) * dt) + s((p - 1) * dt)) + (1 -2 * K) * Z(p - 1, i) + (2 * K) * Z(p - 1, i + 1); 
       elseif i > 1 && i < m + 1
            C(i , i - 1) = K/(4 * (i - 1)) - K/2;
            C(i , i) = 1 + K;
            C(i, i + 1) = -1 * (K/2) -1 * (K /(4 * (i - 1)));
            D(i, 1) = ((a * dt) / (2 * k)) * (s((p - 1 - 1) * dt) + s((p - 1) * dt)) + ((K/2) -1 * (K /(4 * (i - 1)))) * Z(p - 1, i - 1) + (1 - K) * Z(p - 1, i) + ((K/(4 * (i - 1))) + (K/2)) * Z(p - 1, i + 1);
       elseif i == m + 1
            C(m + 1 , m) = -1 * K;
            C(m + 1  , m + 1) = 1 + K + (-1 * K/2 -1 * K/(4*(i-1))) * ((-2 * h * dr)/k);
            D(m + 1, 1) = ((a * dt) / (2 * k)) * (s((p - 1 - 1) * dt) + s((p - 1) * dt)) + (K + K/(2*(i-1))) * ((2 * h * dr * Timbient)/k) + K * Z(p - 1, i -1) + (1 - 1 * K + (K/2 + K/(4*(i-1))) * ((-2 * h * dr)/k)) * Z(p - 1, i);
       end
    end
    B = C \ D;
    Z(p, : ) = B;
end
%% Plotting T vs t
figure('Name', 'Question 2', 'NumberTitle', 'off')
subplot(2, 1, 1)
plot(t, A( : , r1/dr + 1))
grid on
xlabel('Time(s)')
ylabel('Temperature(°C) in r=4mm')
title('T in r=4mm vs t')
hold on
plot(t, Z( : , r1/dr + 1))
grid on
legend('Explicit', 'Crank–Nicolson')
subplot(2, 1, 2)
plot(t, A( : , r2/dr + 1))
grid on
xlabel('Time(s)')
ylabel('Temperature(°C) in r=8mm')
title('T in r=8mm vs t')
hold on
plot(t, Z( : , r2/dr + 1))
grid on
legend('Explicit', 'Crank–Nicolson')
%% Table for selected points
Time = [900, 1800, 3600];
q = [A(t1/dt + 1, r1/dr + 1), A(t2/dt + 1, r1/dr + 1), A(t3/dt + 1, r1/dr + 1)];
r = [Z(t1/dt + 1, r1/dr + 1), Z(t2/dt + 1, r1/dr + 1), Z(t3/dt + 1, r1/dr + 1)];
s = [A(t1/dt + 1, r2/dr + 1), A(t2/dt + 1, r2/dr + 1), A(t3/dt + 1, r2/dr + 1)];
u = [Z(t1/dt + 1, r2/dr + 1), Z(t2/dt + 1, r2/dr + 1), Z(t3/dt + 1, r2/dr + 1)];
T = table(Time', q', r', s', u', 'VariableNames', {'Time', 'r=4mm(Explicit)', 'r=4mm(Crank–Nicolson)', 'r=8mm(Explicit)', 'r=8mm(Crank–Nicolson)'});
disp(T)