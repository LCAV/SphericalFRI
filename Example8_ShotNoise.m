%==========================================================================
% Shot-noise removal (detecting corrupted samples)
%==========================================================================

L       = 6;
L_prime = 12;
K       = 8;

% Generate a random bandlimited signal
f          = rand(2*L, 2*L-1) * 1;
f_lm       = ssht_forward(f, L, 'Method', 'DH', 'Reality', true);

% Need to reproject to really get samples of a bandlimited signal (because
% we don't have critical sampling!)
f          = ssht_inverse(f_lm, L, 'Method', 'DH', 'Reality', true);

% Zero padding
f_prime_lm = zeros(L_prime^2, 1);
f_prime_lm(1:L^2) = f_lm;

f_prime = ssht_inverse(f_prime_lm, L_prime, 'Method', 'DH');

% Generate the sampling grid
[t, p] = ssht_sampling(L_prime, 'Method', 'DH', 'Grid', true);

% Corrupt some samples
idx_corrupt = ceil(rand(K, 1) * numel(f_prime));
s           = rand(K, 1) * 2;
% s = -f_prime(idx_corrupt);
f_corrupt   = f_prime;
f_corrupt(idx_corrupt) = f_corrupt(idx_corrupt) + s;

%--------------------------------------------------------------------------
% Add iid Gaussian noise to measurements
%--------------------------------------------------------------------------

SNR = 200; % SNR in dB

sS2 = var(f_corrupt(:));
sN2 = sS2 / 10^(SNR/10);

noise = randn(size(f_corrupt));
noise = sqrt(sN2) / sqrt(var(noise(:))) * noise;

f_corrupt = f_corrupt + noise;


% Compute the spectrum of the corrupted signal. Note the trap... use the
% matrix method, not the efficient method!!! Argh, you can also use the
% efficient method... just be careful about your definitions of everything.
% f_corrupt_lm = ssht_forward(f_corrupt, L_prime, 'Method', 'DH');

% So first find the matrix

a = WeightsDH(L_prime);
a_grid = a(:) * ones(1, 2*L_prime-1);
A = zeros(L_prime^2, 2*L_prime*(2*L_prime-1));
i = 0;
for l = 0:L_prime-1
    for m = -l:l
        i = i + 1;
        A(i, :) = sqrt(2*pi)/(2*L_prime)*(a_grid(:) .* conj(SphericalHarmonic(l, m, t(:), p(:))))';
    end
end

B = zeros(L_prime^2, 2*L_prime*(2*L_prime-1));
for i = 1:size(B, 2)
    x = zeros(size(t(:)));
    x(i) = 1;
    B(:, i) = ssht_forward(reshape(x, 2*L_prime, 2*L_prime-1), L_prime, 'Method', 'DH', 'Reality', true);
end

% A should be equal to B, figure out why this is not the case... OK, it _is_
% now that I changed the sampling for theta in WeightsDH (but this is
% different than in the paper, it seems to me). Aaah, they have another
% paper...

f_corrupt_lm = B * f_corrupt(:);

% Compute the spectrum that would be generated only by this particular
% /added/ Dirac
f_only = zeros(size(f_corrupt));
f_only(idx_corrupt) = s;
f_only_lm = ssht_forward(f_only, L_prime, 'Method', 'DH', 'Reality', true);
norm([f_corrupt_lm - f_only_lm])
norm([f_corrupt_lm(L^2+1:end) - f_only_lm(L^2+1:end)])

% Ah, these will be the same because both for f_corrupt and for f_only
% computing the transform gives me the projection onto bandlimited
% functions, and the projection is additive, hence equality. Now what I
% need to do is, I really have to compute the low-pass part of the spectrum
% of these Diracs. Sure enough, I should also scale it by a_j^(L_prime).

[j_corrupt, ~] = ind2sub(size(t), idx_corrupt); % This is actually (j+1)
f_only_direct_lm = DiracSpectrum([a(j_corrupt).*s t(idx_corrupt) p(idx_corrupt)], L_prime-1);

% [f_only_lm(L^2+1:end) f_corrupt_lm(L^2+1:end)]

[th_est, phi_est, a_est] = SphereFRI_SN(f_corrupt_lm, K, L);
[sort(th_est) sort(t(idx_corrupt))]
% [sort(phi_est - pi) sort(p(idx_corrupt))]


% Correct the wrong samples using estimated values. Assume that the
% example was run in the correct operating regime, so that the indices of
% corrupted samples are correct.

[~, idxidx] = min(distance(th_est', t(idx_corrupt)'));
real((a_est(idxidx) ./ a(j_corrupt)) ./ s)

const = 3.863376467968997e-01; % HACK. Fix this later.

f_fixed = f_corrupt;
s_est = real((a_est(idxidx) ./ a(j_corrupt)) / const);
f_fixed(idx_corrupt) = f_fixed(idx_corrupt) - s_est;


%% Plotting

L_plot = 150;

f_prime_plot_lm = zeros(L_plot^2, 1);
f_prime_plot_lm(1:length(f_prime_lm)) = f_prime_lm;
f_prime_plot = ssht_inverse(f_prime_plot_lm, L_plot);
error_lm = f_corrupt_lm - f_prime_lm;

error_plot_lm = zeros(L_plot^2, 1);
error_plot_lm(1:length(error_lm)) = error_lm;
error_plot = ssht_inverse(error_plot_lm, L_plot);

r_s   = zeros(3, K);
r_est = zeros(3, K);
r_orig = zeros(3, K);

[r_s(1, :), r_s(2, :), r_s(3, :)] = ssht_s2c(t(idx_corrupt), p(idx_corrupt), 1 + abs(f_corrupt(idx_corrupt)));
[r_est(1, :), r_est(2, :), r_est(3, :)] = ssht_s2c(t(idx_corrupt), p(idx_corrupt), 1 + abs(s_est + f_prime(idx_corrupt)));
[r_orig(1, :), r_orig(2, :), r_orig(3, :)] = ssht_s2c(t(idx_corrupt), p(idx_corrupt), 1 + abs(f_prime(idx_corrupt)));

%--------------------------------------------------------------------------
% Spherical plot of the function without corruptions
%--------------------------------------------------------------------------

figure(1); clf;
hold all;
colormap jet;

ssht_plot_sphere(abs(f_prime_plot), L_plot, 'Type', 'colour', ...
    'ColourBar', true, 'Lighting', true);

for i = 1:K
%     line([0 r_s(1,i)], [0 r_s(2,i)], [0 r_s(3,i)], 'Color', 'b', 'LineWidth', 4);
    line([0 r_orig(1,i)], [0 r_orig(2,i)], [0 r_orig(3,i)], 'Color', 'b', 'LineWidth', 6);
end
scatter3(r_orig(1,:), r_orig(2,:), r_orig(3,:), 200, 'bd', 'filled', 'LineWidth', 2);
axis tight;
cax = caxis;

%--------------------------------------------------------------------------
% Spherical plot of the error function due to corruptions
%--------------------------------------------------------------------------

figure(2); clf;
hold all;

ssht_plot_sphere(abs(error_plot), L_plot, 'Type', 'colour', ...
    'ColourBar', true, 'Lighting', true);

for i = 1:K
    line([0 r_s(1,i)], [0 r_s(2,i)], [0 r_s(3,i)], 'Color', 'b', 'LineWidth', 4);
    line([0 r_orig(1,i)], [0 r_orig(2,i)], [0 r_orig(3,i)], 'Color', 'r', 'LineWidth', 6);
end
scatter3(r_s(1,:), r_s(2,:), r_s(3,:), 200, 'bd', 'filled', 'LineWidth', 2);
axis tight;

caxis(cax);

%--------------------------------------------------------------------------
% Spherical plot of the corrupted function
%--------------------------------------------------------------------------

figure(3); clf;
colormap jet;
hold all;

ssht_plot_sphere(abs(f_prime_plot + error_plot), L_plot, 'Type', 'colour', ...
    'ColourBar', true, 'Lighting', true);

for i = 1:K
    line([0 r_s(1,i)], [0 r_s(2,i)], [0 r_s(3,i)], 'Color', 'b', 'LineWidth', 4);
    line([0 r_est(1,i)], [0 r_est(2,i)], [0 r_est(3,i)], 'Color', 'r', 'LineWidth', 6);
end
scatter3(r_s(1,:), r_s(2,:), r_s(3,:), 200, 'bd', 'filled', 'LineWidth', 2);
scatter3(r_est(1,:), r_est(2,:), r_est(3,:), 200, 'ro', 'LineWidth', 2);
axis tight;

caxis(cax);

%%
%--------------------------------------------------------------------------
% Different 2-D plots
%--------------------------------------------------------------------------

h = figure(4); clf;

col = summer;
col = col(1:25:end, :);

set(h,'NextPlot','replacechildren')
set(h,'DefaultAxesColorOrder', col);

hold all;
plot(abs(f_prime(:)));
stem(idx_corrupt, abs(f_corrupt(idx_corrupt) - f_prime(idx_corrupt)),'d', 'filled', 'MarkerSize', 4, 'LineWidth', 1.2);
stem(idx_corrupt, abs(f_fixed(idx_corrupt) - f_prime(idx_corrupt)), 's', 'filled', 'MarkerSize', 4, 'LineWidth', 1.2);

xlabel('Linearized sample index');
ylabel('Magnitudes');
legend('Bandlimited signals', 'Sample corruption', 'Error after correction');

figWidth = 3;
ExportifyFigure(h, [figWidth, 2/3*figWidth]);

filename = sprintf('../papers/TSP/Figures/ShotNoisePlots.eps');
% export_fig('format','eps', filename);
