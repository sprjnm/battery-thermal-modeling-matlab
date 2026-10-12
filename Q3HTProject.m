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
radius = 0:dr:R;
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
%% Plotting T vs r
figure('Name', 'Question 3 (1)', 'NumberTitle', 'off')
subplot(3, 1, 1)
plot(radius, A(t1/dt + 1 , : ))
grid on
xlabel('Radius(m)')
ylabel('Temperature(°C) in t=900s')
title('T in t=900s vs r')
hold on
plot(radius, Z(t1/dt + 1 , : ))
grid on
legend('Explicit', 'Crank–Nicolson')
subplot(3, 1, 2)
plot(radius, A(t2/dt + 1 , : ))
grid on
xlabel('Radius(m)')
ylabel('Temperature(°C) in t=1800s')
title('T in t=1800s vs r')
hold on
plot(radius, Z(t2/dt + 1 , : ))
grid on
legend('Explicit', 'Crank–Nicolson')
subplot(3, 1, 3)
plot(radius, A(t3/dt + 1 , : ))
grid on
xlabel('Radius(m)')
ylabel('Temperature(°C) in t=3600s')
title('T in t=3600s vs r')
hold on
plot(radius, Z(t3/dt + 1 , : ))
grid on
legend('Explicit', 'Crank–Nicolson')
%% Temperature distribution map
Tt1e = A(t1/dt + 1 , : );
Tt2e = A(t2/dt + 1 , : );
Tt3e = A(t3/dt + 1 , : );
Tt1c = Z(t1/dt + 1 , : );
Tt2c = Z(t2/dt + 1 , : );
Tt3c = Z(t3/dt + 1 , : );
theta = linspace(0, 2*pi, 1000);
[Theta, Radius] = meshgrid(theta, radius);   
[X, Y] = pol2cart(Theta, Radius);
T_2D1e = repmat(Tt1e', 1, length(theta));
T_2D2e = repmat(Tt2e', 1, length(theta));
T_2D3e = repmat(Tt3e', 1, length(theta));
T_2D1c = repmat(Tt1c', 1, length(theta));
T_2D2c = repmat(Tt2c', 1, length(theta));
T_2D3c = repmat(Tt3c', 1, length(theta));
figure('Name', 'Question 3 (2)', 'NumberTitle', 'off')
subplot(3,2,1)
pcolor(X * 1000, Y * 1000, T_2D1e);  
shading interp;                  
colormap(jet);                    
colorbar;                        
title(['Temperature distribution(°C) at t = ', num2str(t1), 's (Explicit)']);
xlabel('X (mm)');
ylabel('Y (mm)');
axis equal;
subplot(3,2,3)
pcolor(X * 1000, Y * 1000, T_2D2e);  
shading interp;                  
colormap(jet);                    
colorbar;                        
title(['Temperature distribution(°C) at t = ', num2str(t2), 's (Explicit)']);
xlabel('X (mm)');
ylabel('Y (mm)');
axis equal;
subplot(3,2,5)
pcolor(X * 1000, Y * 1000, T_2D3e);  
shading interp;                  
colormap(jet);                    
colorbar;                        
title(['Temperature distribution(°C) at t = ', num2str(t3), 's (Explicit)']);
xlabel('X (mm)');
ylabel('Y (mm)');
axis equal;
subplot(3,2,2)
pcolor(X * 1000, Y * 1000, T_2D1c);  
shading interp;                  
colormap(jet);                    
colorbar;                        
title(['Temperature distribution(°C) at t = ', num2str(t1), 's (Crank–Nicolson)']);
xlabel('X (mm)');
ylabel('Y (mm)');
axis equal;
subplot(3,2,4)
pcolor(X * 1000, Y * 1000, T_2D2c);  
shading interp;                  
colormap(jet);                    
colorbar;                        
title(['Temperature distribution(°C) at t = ', num2str(t2), 's (Crank–Nicolson)']);
xlabel('X (mm)');
ylabel('Y (mm)');
axis equal;
subplot(3,2,6)
pcolor(X * 1000, Y * 1000, T_2D3c);  
shading interp;                  
colormap(jet);                    
colorbar;                        
title(['Temperature distribution(°C) at t = ', num2str(t3), 's (Crank–Nicolson)']);
xlabel('X (mm)');
ylabel('Y (mm)');
axis equal;