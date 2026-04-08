E0 = 1.0;
t0 = 10;
t_center = 10; 
c = 1.0;
lambda0 = 2;
k = 2 * pi / lambda0;

mode = 4;

if mode == 1
    n = 1.0;
    dn_dlambda = -0.1;
    w = 2*pi*c/lambda0;
    vp = w/k;
    vg = c/(n-lambda0*dn_dlambda);
elseif mode == 2
    w = 2*pi*c/lambda0;
    vp = w/k;
    vg = 0;
elseif mode == 3
    w = 0;
    vp = 0;
    vg = 2.0;
elseif mode == 4
    n = 1.5;
    dn_dlambda = 0.2;
    w = 2*pi*c/lambda0;
    vg = c/(n-lambda0*dn_dlambda);
    vp = -vg;

end

fprintf('Phase velocity (vp) = w/k = %.2f\n', vp);
fprintf('Group velocity (vg) = c/n - lambda0*dn_dlambda = %.2f\n', vg);

z = linspace(0, 100, 500); 
t_vector = linspace(0,100,100);  

fig = figure('Name', 'Gaussian pulse', 'Position', [100, 100, 950, 600]);
ax = axes(fig);
hold on; grid on;

b_wave = plot(z, zeros(size(z)), 'b', 'LineWidth', 1.2, 'DisplayName', 'Electric Field');
r_env = plot(z, zeros(size(z)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Envelope');

point_group = plot(0, 0, 'ro', 'Markersize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'Group point (vg)');
point_phase = plot(0, 0, 'go', 'Markersize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 1, 'DisplayName', 'Phase point (vp)');

axis([min(z), max(z), -1.3* E0, 1.3* E0]); 
xlabel('Position (t)');
ylabel('Electric Field (E)');
legend('show', 'Location', 'northeast');

while ishandle(fig) 
    for i = 1:length(t_vector)
    t1 = t_vector(i); 
    
    t_env = z - t_center - vg * t1; 
    phi_carrier = k * (z - t_center) - w * t1; 
    
    envelope = E0 * exp(-t_env.^2 / (2 * t0^2));
    E = envelope .* cos(phi_carrier);
    
    t_group_point = t_center + vg * t1;
    y_group_point = E0;
    
    t_phase_point = t_center + vp * t1;
    y_phase_point = E0 * exp(-(t_phase_point - t_center - vg * t1).^2 / (2 * t0^2));
    
    set(r_env, 'YData', envelope); 
    set(b_wave, 'YData', E);
    set(point_group, 'XData', t_group_point, 'YData', y_group_point);
    set(point_phase, 'XData', t_phase_point, 'YData', y_phase_point);
    
    drawnow;
    pause(0.02);
    end

   if ishandle(fig)
    pause(0.5);
   end
end









