%==========================================================================
% Example of diffusion reconstruction
%==========================================================================

L_cutoff = 7;     % Cutoff bandwidth for spectral computation (0...L_cutoff-1)
L_plot   = 100;   % Bandwidth for plotting (to get nice figures)
k        = 0.01;   % Diffusion constant
t0       = 1;     % Time to sample
K        = 3;     % Number of sources

% Sampling grid
[t, p] = ssht_sampling(L_cutoff, 'Grid', true);

% Compute the diffusion kernel in spectral and spatial domain
l = (0:L_cutoff-1)';
u_lm = zeros(L_cutoff^2, 1);
u_lm(l.^2 + l + 1) = exp(-l.*(l+1)*k*t0);
u = ssht_inverse(u_lm, L_cutoff, 'Reality', true);

% Generate random source locations
mag_s = rand(K, 1)*3 + 1;
phi_s = rand(K, 1)*2*pi;
th_s  = pi/2 - asin((rand(K, 1)*2-1));

% Compute the spectrum of Diracs
fprintf('Computing the spectrum of Diracs... ');
f_lm = sum(DiracSpectrum([mag_s th_s phi_s], L_cutoff), 2);
fprintf('Done!\n');

% Do the convolution
l_lin = round(real(sqrt(4*(1:L_cutoff^2) - 3)/2 - 1/2))'; % A trick (document this somewhere)
g_lm  = sqrt(1./(2*l_lin+1)) .* exp(-l_lin.*(l_lin+1)*k*t0) .* f_lm;

% Compute the function in the spatial domain
g = ssht_inverse(g_lm, L_cutoff, 'Reality', true);

% Add iid Gaussian noise to measurements

SNR = 50; % SNR in dB

sS2 = var(g(:));
sN2 = sS2 / 10^(SNR/10);
noise = randn(size(g));
noise = sqrt(sN2) / sqrt(var(noise(:))) * noise;
g = g + noise;
g_lm = ssht_forward(g, L_cutoff, 'Reality', true);

% Perform the deconvolution by spectral division (use the \ell trick again
% here, and also in other example from where you copied this)
LK = ceil(K + sqrt(K));
f_hat = zeros(LK^2, 1);
i = 0;
for ell = 0:LK-1
    for m = -ell:ell
        i = i + 1;
        f_hat(i) = sqrt(2*ell+1)*g_lm(i) / u_lm(ell^2 + ell + 1);
    end
end

[th_est, phi_est, mag_est] = SphereFRI(f_hat, K);
mag_est = real(mag_est);

%% Plotting

close all;

%--------------------------------------------------------------------------
% Sphere with real and estimated peaks
%--------------------------------------------------------------------------

figure(1);
hold on;

g_plot_lm = zeros(L_plot^2, 1);
g_plot_lm(1:length(g_lm)) = g_lm;
g_plot = ssht_inverse(g_plot_lm, L_plot, 'Reality', true);

ssht_plot_sphere(g_plot, L_plot, 'Type', 'colour', ...
    'ColourBar', false, 'Lighting', true);

% material dull;
% alpha(0.7);

r_s   = zeros(3, K);
r_est = zeros(3, K);

[r_s(1, :), r_s(2, :), r_s(3, :)] = SphCart(th_s, phi_s, mag_s);
[r_est(1, :), r_est(2, :), r_est(3, :)] = SphCart(th_est, phi_est, mag_est);

for i = 1:K
    line([0 r_s(1,i)], [0 r_s(2,i)], [0 r_s(3,i)], 'LineWidth', 2);
    line([0 r_est(1,i)], [0 r_est(2,i)], [0 r_est(3,i)], 'Color', 'r', 'LineWidth', 3);
end

scatter3(r_s(1,:), r_s(2,:), r_s(3,:), 50, 'bd', 'filled', 'LineWidth', 2);
scatter3(r_est(1,:), r_est(2,:), r_est(3,:), 100, 'ro', 'LineWidth', 2);

axis([-4 4 -4 4 -4 4])
colormap jet;

