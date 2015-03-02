%==========================================================================
%  This example plots Helmholtz Green's functions on the sphere
%==========================================================================

% NOTE: It would be more efficient to directly compute the spectral
% coefficients, but this way we are trying to emulate what's happpening in
% the real world.

% Some parameters
L = 100;            % Cutoff bandwidth for spectral computation (0...L-1)
Lcutoff = 30;       % Cutoff bandwidth for Green's function computation (0...L-1)
r = 0.25;           % Rigid sphere radius
f = 500;            % Frequency in Hertz
k = 2*pi*f/343;     % Wavenumber (== 2*pi*f/c)
r_s = [0 0 1;
       0 0 2; 
       0 0 3; 
       0 0 4; 
       0 0 5]';

% Sampling grid
[t, p] = ssht_sampling(L, 'Grid', true);
[x, y, z] = SphCart(t, p, r*ones(size(t)));

g = zeros([size(t, 1) size(t, 2) size(r_s, 2)]);
g_lm = zeros(L^2, size(r_s, 2));

for i = 1:size(r_s, 2)
    % Compute the Green's function
    disp('Computing Green''s Function...');
    g_i = GreenRigidSphere([x(:) y(:) z(:)]', r_s(:, i), k, Lcutoff);
    g(:, :, i) = reshape(g_i, size(t));P

    % And its Spherical Harmonic Transform
    disp('Computing SHT...');
    g_lm(:, i) = ssht_forward(g(:, :, i), L);
end


%% Do the plots

% Green's function on the sphere
figure(1);
subplot(2,2,1);
ssht_plot_sphere(abs(g(:, :, 1)), L, 'Type', 'colour', ...
    'ColourBar', true, 'Lighting', true);

subplot(2,2,2);
ssht_plot_sphere(abs(g(:, :, 2)), L, 'Type', 'colour', ...
    'ColourBar', true, 'Lighting', true);

subplot(2,2,3);
ssht_plot_sphere(abs(g(:, :, 3)), L, 'Type', 'colour', ...
    'ColourBar', true, 'Lighting', true);

subplot(2,2,4);
ssht_plot_sphere(abs(g(:, :, 4)), L, 'Type', 'colour', ...
    'ColourBar', true, 'Lighting', true);


% Green's function spectrum
Lplot = 20;
l = 0:Lplot;
l0 = l.^2 + l + 1;

figure(2);
clf;
hold all;

stem(0:length(l0)-1, abs(g_lm(l0, 1)), 'bv', 'filled', 'MarkerSize', 6, 'LineWidth', 1.2);
stem(0:length(l0)-1, abs(g_lm(l0, 2)), 'rs', 'filled', 'MarkerSize', 6, 'LineWidth', 1.2);
stem(0:length(l0)-1, abs(g_lm(l0, 3)), 'ko', 'filled', 'MarkerSize', 6, 'LineWidth', 1.2);
stem(0:length(l0)-1, abs(g_lm(l0, 4)), 'mp', 'filled', 'MarkerSize', 6, 'LineWidth', 1.2);
ax = axis;
grid;
xlabel('l');
ylabel('|G(l, 0)|');

g_aux = g(:, p==0, 1:4);
g_aux = squeeze(g_aux(:, 1, :));
pbaspect([1 1 1]);

% Green's function along a meridian
s = figure(3);
clf;
hold all;

t0 = t(p==0);
plot(t0, abs(g_aux(:, 1)), 'b--', 'LineWidth', 1.2);
plot(t0, abs(g_aux(:, 2)), 'r-.', 'LineWidth', 1.2);
plot(t0, abs(g_aux(:, 3)), 'k:',  'LineWidth', 1.2);
plot(t0, abs(g_aux(:, 4)), 'm-',  'LineWidth', 1.2);

sp = 10;
plot(t0(1:sp:end), abs(g_aux(1:sp:end, 1)), 'bv', 'LineWidth', 1.5);
plot(t0(1:sp:end), abs(g_aux(1:sp:end, 2)), 'rs', 'LineWidth', 1.5);
plot(t0(1:sp:end), abs(g_aux(1:sp:end, 3)), 'ko',  'LineWidth', 1.5);
plot(t0(1:sp:end), abs(g_aux(1:sp:end, 4)), 'mp',  'LineWidth', 1.5);

legend('1 m', '2 m', '3 m', '4 m');

axis tight;
pbaspect([1 .7 1]);
grid;
xlabel('\theta');
ylabel('g');
