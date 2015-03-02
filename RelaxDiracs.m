function [th_r, ph_r, mag_r] = RelaxDiracs(th, ph, mag, f_lm)

L = sqrt(numel(f_lm));


S0 = [mag(:) th(:) ph(:)];

opt = optimset('MaxFunEval', 1e6, ...
               'MaxIter',    1e6);
S = fminsearch(@(S) sum( abs(DiracSpectrum(S, L) - f_lm).^2 ), S0, opt);

th_r = S(:, 2);
ph_r = S(:, 3);
mag_r = S(:, 1);


