clear all; close all; clc;

T0 = 1e-12;          
t1 = 10*T0;       
Nt = 1024;          
t = linspace(-t1, t1, Nt);
dt = t(2) - t(1);

z2 = 10;         
Nz = 100;           
z = linspace(0, z2, Nz);

beta2 = -20e-27;    
LD = T0^2/abs(beta2); 
z1 = z * LD;     
P0 = 1;              
C = 0;               

p1 = 1; 
m = 3; 

switch p1
    case 1
        A0 = sqrt(P0) * exp(-0.5*(t/T0).^2);
    case 2
        A0 = sqrt(P0) * exp(-0.5*(t/T0).^(2*m));
    case 3
        A0 = sqrt(P0) * sech(t/T0);
    case 4
        A0 = sqrt(P0) * double(abs(t) <= T0);
end

omega = 2*pi*(-Nt/2:Nt/2-1)/(Nt*dt);
omega = fftshift(omega);
dz = z(2) - z(1);

A = zeros(Nz, Nt);
A(1,:) = A0;

if p1 == 1
    for k = 2:Nz
        zeta = z(k);
        factor = 1 + 1i*sign(beta2)*zeta;
        A(k,:) = sqrt(P0) / sqrt(factor) * exp(-0.5*(1+1i*C)*(t/T0).^2 ./ factor);
    end
else
    A(2:end, :) = NaN; 
end

A1 = zeros(Nz, Nt);
A1(1,:) = A0;

for k = 2:Nz
    D1 = exp(1i * 0.5 * beta2 * omega.^2 * dz * LD);
    A2 = fft(A1(k-1,:));
    A2 = A2 .* D1;
    A1(k,:) = ifft(A2);
end

T1 = T0 * sqrt(1 + (z*sign(beta2)).^2);

figure('Position', [100, 100, 1200, 800]);

subplot(2,3,1);
plot(t/T0, abs(A0).^2, 'b-', 'LineWidth', 2);
xlabel('Normalized time t/T_0');
ylabel('Power |A|^2');
title('Initial Pulse Shape');
grid on;

subplot(2,3,2);
imagesc(t/T0, z, abs(A).^2);
xlabel('Normalized time t/T_0');
ylabel('Normalized distance z/L_D');
title('Pulse Evolution (Analytical)');
colorbar;
colormap('jet');

subplot(2,3,3);
z3 = [1, round(Nz/3), round(2*Nz/3), Nz];
hold on;
colors = ['b', 'r', 'g', 'm'];
for k = 1:length(z3)
    idx = z3(k);
    plot(t/T0, abs(A1(idx,:)).^2, colors(k), 'LineWidth', 1.5);
end
xlabel('Normalized time t/T_0');
ylabel('Power |A|^2');
title('Pulse at Different Distances (Num)');
legend('z/LD=0', 'z/LD=3.3', 'z/LD=6.7', 'z/LD=10');
grid on;

subplot(2,3,4);
if p1 == 1
    plot(z, T1/T0, 'b-', 'LineWidth', 2);
    title('Pulse Broadening (Theoretical)');
else
    title('Broadening Theory N/A for this pulse');
end
xlabel('Normalized distance z/L_D');
ylabel('Normalized pulse width T(z)/T_0');
grid on;

subplot(2,3,5);
z4 = Nz; 
plot(t/T0, abs(A(z4,:)).^2, 'b-', 'LineWidth', 2);
hold on;
plot(t/T0, abs(A1(z4,:)).^2, 'r--', 'LineWidth', 2);
xlabel('Normalized time t/T_0');
ylabel('Power |A|^2');
title(['Comparison at z = ', num2str(z(z4)), ' L_D']);
legend('Analytical', 'Numerical (FFT)');
grid on;

subplot(2,3,6);
imagesc(t/T0, z, abs(A1).^2);
xlabel('Normalized time t/T_0');
ylabel('Normalized distance z/L_D');
title('Pulse Evolution (Numerical)');
colorbar;
colormap('jet');