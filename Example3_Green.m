%==========================================================================
%  This example plots parts of Fig. 7 from the paper (C, F == the source
%  localization example)
%==========================================================================


% Some parameters
L = 100;            % Cutoff bandwidth for spectral computation (0...L-1)
Lcutoff = 30;       % Cutoff bandwidth for Green's function computation (0...Lcutoff-1)
r = 0.20;           % Rigid sphere radius
k = 2*pi*3000/343;  % Wavenumber (== 2*pi*f/c)
K = 3;

% Sampling grid
[t, p] = ssht_sampling(L, 'Grid', true);
[x, y, z] = SphCart(t, p, r*ones(size(t)));


% Compute the filter spectrum
h = GreenRigidSphere([x(:) y(:) z(:)]', [0 0 4]', k, Lcutoff);
h = reshape(h, size(t));
h_lm = ssht_forward(h, L);

% Generate random source locations
rad_s = rand(K, 1)*2 + 3;
phi_s = rand(K, 1)*2*pi;
th_s  = pi/2 - asin((rand(K, 1)*2-1));

r_s   = zeros(3, K);
[r_s(1, :), r_s(2, :), r_s(3, :)] = SphCart(th_s, phi_s, rad_s);

% Generate the spectrum for every angular spike. We deliberately don't use
% convolution here, because we want to take into account the variation of
% the Green's function with the source's distance.
g = zeros([size(t, 1) size(t, 2) K]);
g_lm = zeros(L^2, K);
for i = 1:K
    % Compute the Green's function
    disp('Computing Green''s function...');
    g_i = GreenRigidSphere([x(:) y(:) z(:)]', r_s(:, i), k, Lcutoff);
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

% GreenRigidSphere uses a different coordinate system, so compensate
% Figure out for sure where the pi comes from (some multiplication by (-1)^m)
% Check out also Example7 for this.

phi_est = mod(phi_est, 2*pi);

[sort(th_s) sort(th_est)]

%% Do the plots
figure(1); clf;
ssht_plot_sphere(abs(sum(g, 3)), L, 'Type', 'colour', 'Lighting', true);
hold all;

scale = 1;
r_s = r_s * scale;

scatter3(r_s(1,:), r_s(2,:), r_s(3,:), 300, 'd', 'filled', 'LineWidth', 3);

r_est = zeros(size(r_s));

[r_est(1, :), r_est(2, :), r_est(3, :)] = SphCart(th_est, phi_est, 2);
for k = 1:K
    line([0 r_s(1,k)], [0 r_s(2,k)], [0 r_s(3,k)], 'LineWidth', 1);
    line([0 r_est(1,k)], [0 r_est(2,k)], [0 r_est(3,k)], 'Color', 'r', 'LineWidth', 4);
end

axis tight equal;

figure(2); clf;
ssht_plot_sphere(abs(sum(g, 3)), L, 'Type', 'colour', ...
    'ColourBar', true, 'Lighting', true, 'type', 'parametric');

