%--------------------------------------------------------------------------
% CRLB example for 1 Dirac
%--------------------------------------------------------------------------

% Dirac parameters
t0_list = [pi/2 pi/4];
p0_list = [1 1];
a0_list  = [1 0.1];

% Sampling grid
L = 2;
[T, P] = ssht_sampling(L, 'Grid', true);

% Simulation parameters
SNR_range = -10:1:50;
n_iter    = 100;


% Compute the samples of the bandlimited Dirac

% Simulate and plot
for ifig = 1:2
    h = figure(ifig);
    col = summer;
    col = col(1:32:end, :);
    set(h,'NextPlot','replacechildren')
    set(h,'DefaultAxesColorOrder', col);
    clf; hold all;
end

matlabpool(3);

var_out = zeros(length(t0_list), length(SNR_range));
var_theta = length(t0_list), length(SNR_range);
var_noise = length(t0_list), length(SNR_range);
CRLB = length(t0_list), length(SNR_range);
CRLB_theta = length(t0_list), length(SNR_range);

for i_pos = 1:length(t0_list)
    t0 = t0_list(i_pos);
    p0 = p0_list(i_pos);
    a0 = a0_list(i_pos);
        
    flm = DiracSpectrum([a0, t0, p0], L);
    f   = real(ssht_inverse(flm, L));

    % Signal power
    P_signal = var(f(:));
    
  
    parfor i_SNR = 1:length(SNR_range)
        SNR = SNR_range(i_SNR)
        estimates = zeros(n_iter, 3);
        for iter = 1:n_iter
            % Get the samples of the bandlimited Dirac
            P_noise = P_signal / 10^(SNR/10);
            noise   = randn(size(f));
            noise   = sqrt(P_noise / var(noise(:))) * noise;

            f_noise = f + noise;
            flm_noise = ssht_forward(f_noise, L, 'Reality', true);
            [t_est, p_est, a_est] = SphereFRI(flm_noise, 1, false, true);
            estimates(iter, :) = [t_est, p_est, a_est];
        end
        var_out(i_pos, i_SNR) = sum(sum(bsxfun(@minus, real(estimates(:, 1:2)), [t0 p0]).^2))/n_iter; 
        var_theta(i_pos, i_SNR) = sum(real(estimates(:, 1) - t0).^2)/n_iter; 
        var_noise(i_pos, i_SNR) = P_noise;
        
        F_inv = CramerRaoBound(T(:), P(:), P_noise, t0, p0, a0, L, 1);
        CRLB(i_pos, i_SNR) = (F_inv(1,1) + F_inv(2,2));
        CRLB_theta(i_pos, i_SNR) = F_inv(1,1);
    end

    % Here you do the variance of the estimator... You should instead do the
    % MSE. They interestingly coincide... so unbiased...
    
    % Figure 1
    figure(1);
    plot(SNR_range, 10*log10(var_out(i_pos, :)), ... 
        ':', ...
        'LineWidth', 1, ...
        'MarkerSize', 3, ...
        'MarkerFaceColor', [0.5,0.5,0.5], ...
        'MarkerEdgeColor', 'k');
    plot(SNR_range, 10*log10(CRLB(i_pos, :)), ...
            '-', ...
            'LineWidth', 1, ...
            'MarkerSize', 3, ...
            'MarkerFaceColor', [0.5,0.5,0.5], ...
            'MarkerEdgeColor', 'k');
        
    % Figure 2 (only theta)
    figure(2);
    plot(SNR_range, 10*log10(var_theta(i_pos, :)), ... 
        ':', ...
        'LineWidth', 1, ...
        'MarkerSize', 3, ...
        'MarkerFaceColor', [0.5,0.5,0.5], ...
        'MarkerEdgeColor', 'k');
    plot(SNR_range, 10*log10(CRLB_theta(i_pos, :)), ...
            '-', ...
            'LineWidth', 1, ...
            'MarkerSize', 3, ...
            'MarkerFaceColor', [0.5,0.5,0.5], ...
            'MarkerEdgeColor', 'k');
    filename = sprintf('CRLB_data_%d.mat', i_pos);
    save(filename);
end

xlabel('Input SNR [dB]');
ylabel('10log10(Output MSE)');
legend('MSE, \pi/2', 'CRLB, \pi/2', 'MSE, \pi/4', 'CRLB, \pi/4');
axis tight;

figWidth = 3.5;
ExportifyFigure(h, [figWidth, 1/2*figWidth]);

filename = sprintf('../papers/TSP/Figures/CRLB_plots.pdf');
% export_fig('format','eps', filename);

