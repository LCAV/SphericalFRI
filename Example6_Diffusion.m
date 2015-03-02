%==========================================================================
% This example plots the diffusion kernel, and aliasing error
%==========================================================================

L_cutoff = 100; % Cutoff bandwidth for kernel computation
K        = 3;   % Number of Diracs
[t, p]   = ssht_sampling(L_cutoff, 'Grid', true); % Sampling grid


k     =  10.^([-2 -1.5 -1 0]); % Diffusion constant
t0    =  1;               % Time to sample

% Compute the diffusion kernel in spectral and spatial domains
l = (0:L_cutoff-1)';
u_lm = zeros(L_cutoff^2, numel(k));
u_lm(l.^2 + l + 1, :) = exp(-l.*(l+1)*k*t0);

u = zeros([size(t,1) size(t,2) numel(k)]);
for i = 1:numel(k)
    u(:, :, i) = ssht_inverse(u_lm(:, i), L_cutoff, 'Reality', true);
end    
    
% Aliasing error computation
L_alias = 15;
aliasing = zeros(L_alias, numel(k));
for L = 0:L_alias-1
    aliasing(L+1, :) = sum(diag(1./(2*l(L+1:end)+1))*exp(-l(L+1:end).*(l(L+1:end)+1)*k*t0).^2) ...
                    ./ sum(diag(1./(2*l(1:end)+1))*exp(-l(1:end).*(l(1:end)+1)*k*t0).^2);
end


%% Do the plots

close all;

%--------------------------------------------------------------------------
% Relative aliasing error
%--------------------------------------------------------------------------
h = figure(1);

col = autumn;
col = col(1:20:end, :);

% col = [0 0 0;
%        1 0 0;
%        0 0 1];

set(h,'NextPlot','replacechildren')
set(h,'DefaultAxesColorOrder', col);

plot(0:L_alias-1, 10*log10(aliasing), ...
        '-s', ...
        'LineWidth', .4, ...
        'MarkerSize', 3, ...
        'MarkerFaceColor', [0.5,0.5,0.5], ...
        'MarkerEdgeColor', 'k');

xlabel('Cutoff degree l');
ylabel('log(Relative aliasing energy)');
axis tight;

% legend('k = 0.01', 'k = 0.03', 'k = 0.10', 'Location', 'SouthWest');


figWidth = 3;
ExportifyFigure(h, [figWidth, 2/3*figWidth]);

filename = sprintf('../papers/TSP/Figures/RelativeAliasingError.eps');
export_fig('format','eps', filename);

%--------------------------------------------------------------------------
% Shape of the diffusion kernel
%--------------------------------------------------------------------------
h = figure(2);

set(h,'NextPlot','replacechildren')
set(h,'DefaultAxesColorOrder', col);

hold all;

u_aux = u(:, p==0, :);
u_aux = squeeze(u_aux(:, 1, :));

th0 = t(p==0);

for i = 1:numel(k)
    plot(th0, u_aux(1:end-1, i) / u_aux(1, i), 'LineWidth', 1);
end

axis tight;

xlabel('Colatitude \theta');
ylabel('Normalized diffusion kernel');

% legend('k = 0.01', 'k = 0.03', 'k = 0.10', 'Location', 'NorthEast');

figWidth = 3;
ExportifyFigure(h, [figWidth, 2/3*figWidth]);

filename = sprintf('../papers/TSP/Figures/DiffusionKernels.eps');
% export_fig('format','eps', filename);

