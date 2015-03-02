%==========================================================================
%  This example plots parts of Figure 8 from the paper
%==========================================================================

%% -------------------------------------------------------------------------
% Real and imaginary part of the Green's function in spatial and spectral
% domain for two different frequencies.
%--------------------------------------------------------------------------

% Some parameters
L = 100;            % Cutoff bandwidth for spectral computation (0...L-1)
Lcutoff = 30;       % Cutoff bandwidth for Green's function computation (0...L-1)
r = 0.20;           % Rigid sphere radius
freq = [1e3 4e3];   % Frequencies in Hertz
k = 2*pi*freq/343;  % Wavenumber (== 2*pi*f/c)
r_s = [0 0 3]';     % Source position

% Sampling grid
[t, p] = ssht_sampling(L, 'Grid', true);
[x, y, z] = SphCart(t, p, r*ones(size(t)));

g = zeros([size(t, 1) size(t, 2) size(r_s, 2)]);
g_lm = zeros(L^2, size(r_s, 2));

for i = 1:numel(k)
    % Compute the Green's function
    disp('Computing Green''s Function...');
    g_i = GreenRigidSphere([x(:) y(:) z(:)]', r_s, k(i), Lcutoff);
    g(:, :, i) = reshape(g_i, size(t));

    % And its Spherical Harmonic Transform
    disp('Computing SHT...');
    g_lm(:, i) = ssht_forward(g(:, :, i), L);
end

%% ------------------------------------------------------------------------
% Plotting the Green's functions
%--------------------------------------------------------------------------

figWidth = 3;
for i = 1:numel(k)
    % Green's function spectrum
    Lplot = 25;
    l     = 0:Lplot;
    l0    = l.^2 + l + 1;

    h = figure(1);
    clf;
    hold all;

    stem(0:length(l0)-1, real(g_lm(l0, i)), 'bv', 'filled', 'MarkerSize', 2, 'LineWidth', 1);
    stem(0:length(l0)-1, imag(g_lm(l0, i)), 'rs', 'filled', 'MarkerSize', 2, 'LineWidth', 1);
    stem(0:length(l0)-1, abs(g_lm(l0,  i)), 'ko', 'filled', 'MarkerSize', 2, 'LineWidth', 1);
    
    legend('Real', 'Imag', 'Abs');
    if i == 2
        legend('Location', 'SouthEast');
    end
    
    
    axis tight;
    xlabel('Spherical harmonic degree');
    ylabel('|g_l^0|');
    
    % Green's function along a meridian
    g_aux = g(:, p==0, i);

    h = figure(2);
    clf;
    hold all;

    t0 = t(p==0);
    plot(t0, real(g_aux(:, i)), 'b-', 'LineWidth', 0.5);
    plot(t0, imag(g_aux(:, i)), 'r-', 'LineWidth', 0.5);
    plot(t0,  abs(g_aux(:, i)), 'k-',  'LineWidth', 2);

    sp = 5;

    legend('Real', 'Imag', 'Abs', 'Location', 'SouthEast');

    ax = axis;
    ax(2) = pi;
    axis(ax);
    xlabel('Colatitude \theta');
    ylabel('g(\theta, 0)');
end

%% -------------------------------------------------------------------------
% Reconstruction examples at different frequencies
%--------------------------------------------------------------------------

% Some parameters
L  = 12;             % Cutoff bandwidth for spectral computation (0...L-1)
Lcutoff = 30;        % Cutoff bandwidth for Green's function computation (0...L)
r  = 0.20;           % Rigid sphere radius
f0 = 1000;           % Frequency
k0 = 2*pi*f0/343;    % Wavenumber (== 2*pi*f/c)
K  = 3;              % Number of sources to localize (at different frequencies)

% Sampling grid
[t, p] = ssht_sampling(L, 'Grid', true);
[x, y, z] = SphCart(t, p, r*ones(size(t)));

% Compute the filter spectrum
h = GreenRigidSphere([x(:) y(:) z(:)]', [0 0 4]', k0, Lcutoff);
h = reshape(h, size(t));
h_lm = ssht_forward(h, L);

% Generate random source locations
rad_s = rand(K, 1)*2 + 3;
phi_s = rand(K, 1)*2*pi;
th_s  = pi/2 - asin((rand(K, 1)*2-1));

r_s   = zeros(3, K);
[r_s(1, :), r_s(2, :), r_s(3, :)] = SphCart(th_s, phi_s, rad_s);

% Generate the spectrum for every angular spike
g = zeros([size(t, 1) size(t, 2) K]);
g_lm = zeros(L^2, K);
for i = 1:K
    % Compute the Green's function
    disp('Computing Green''s function...');
    g_i = GreenRigidSphere([x(:) y(:) z(:)]', r_s(:, i), k0, Lcutoff);
    g(:, :, i) = reshape(g_i, size(t));

    % And its Spherical Harmonic Transform
    disp('Computing SHT...');
    g_lm(:, i) = ssht_forward(g(:, :, i), L);
end

g_lm_total = sum(g_lm, 2);

% Perform the deconvolution by spectral division
LK = ceil(K + sqrt(K) - 1);
f = zeros((LK+1)^2, 1);
i = 0;
for l = 0:LK
    for m = -l:l
        i = i + 1;
        f(i) = sqrt(2*l+1)*g_lm_total(i) / h_lm(l^2 + l + 1);
    end
end

[th_est, phi_est, a_est] = SphereFRI(f, K);
phi_est = mod(phi_est, 2*pi);

%% ------------------------------------------------------------------------
% Visualizing localization on the sphere
%--------------------------------------------------------------------------
h = figure(3); clf;
ssht_plot_sphere(abs(sum(g, 3)), L, 'Type', 'colour', 'Lighting', true);
hold all;

scale = 1;
r_s = r_s * scale;

scatter3(r_s(1,:), r_s(2,:), r_s(3,:), 50, 'd', 'filled', 'LineWidth', 1);

r_est = zeros(size(r_s));

[r_est(1, :), r_est(2, :), r_est(3, :)] = SphCart(th_est, phi_est, 2);
for k = 1:K
    line([0 r_s(1,k)], [0 r_s(2,k)], [0 r_s(3,k)], 'LineWidth', 1);
    line([0 r_est(1,k)], [0 r_est(2,k)], [0 r_est(3,k)], 'Color', 'r', 'LineWidth', 4);
end

axis tight equal;
