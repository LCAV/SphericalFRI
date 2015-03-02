%==========================================================================
%  This example plots parts of Figure ? from the paper
%==========================================================================

% Some parameters
L = 100;            % Cutoff bandwidth for spectral computation (0...L-1)
Lcutoff = 30;       % Cutoff bandwidth for Green's function computation (0...L-1)
r = 0.20;           % Rigid sphere radius
f = 700;            % Frequency in Hertz
k = 2*pi*f/343;     % Wavenumber (== 2*pi*f/c)
r_s = [0 0 1.0;
       0 0 1.2; 
       0 0 1.4; 
       0 0 1.6; 
       0 0 1.8;
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
    g(:, :, i) = reshape(g_i, size(t));

    % And its Spherical Harmonic Transform
    disp('Computing SHT...');
    g_lm(:, i) = ssht_forward(g(:, :, i), L);
end

%% ------------------------------------------------------------------------
% Plotting the Green's function, and the ratio of Green's functions
%--------------------------------------------------------------------------

% Green's function on the sphere
h = figure(1);
ssht_plot_sphere(abs(g(:, :, 1)), L, 'Type', 'colour', ...
    'ColourBar', false, 'Lighting', true);

figWidth = 2;
ExportifyFigure(h, [figWidth, 3/4*figWidth]);
filename = sprintf('../papers/TSP/Figures/SingleSourceSphere-%d.png', f);
export_fig('format','png', 'renderer', 'opengl', '-m5', filename);


% Spectral ratios
Lplot = 10;
l = 0:Lplot;
l0 = l.^2 + l + 1;

h = figure(2);
clf;
hold all;

for i = 1:size(r_s, 2);
    plot(0:length(l0)-1, abs(g_lm(l0, i)) ./ abs(g_lm(l0, end)), ...
        'k-s', ...
        'LineWidth', .4, ...
        'MarkerSize', 3, ...
        'MarkerFaceColor', [0.5,0.5,0.5]);
end
ax = axis;
grid off;

xlabel('Spherical harmonic degree');
ylabel('Ratio');

figWidth = 2;
ExportifyFigure(h, [figWidth, 2/3*figWidth]);
filename = sprintf('../papers/TSP/Figures/GreenSpectralRatios-%d.eps', f);
export_fig('format','eps', filename);

% Spatial ratios
h = figure(3);
clf;
hold all;

g_aux = g(:, p==0, :);
g_aux = squeeze(g_aux(:, 1, :));
pbaspect([1 1 1]);

t0 = t(p==0);

for i = 1:size(r_s, 2);
    plot(t0, abs(g_aux(:, i)) ./ abs(g_aux(:, end)), 'k-', 'LineWidth', 1);
end

xlabel('Colatitude \theta');
ylabel('Ratio');
axis tight;

ExportifyFigure(h, [figWidth, 2/3*figWidth]);
filename = sprintf('../papers/TSP/Figures/GreenSpatialRatios-%d.eps', f);
export_fig('format','eps', filename);
