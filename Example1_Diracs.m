% Setup some parameteres
K = 5;
L = ceil(K + sqrt(K));

% Generate a random stream of Diracs, and compute its spectrum
S = RandomDiracs(K, 1);

f = DiracSpectrum(S, L);

% Reconstruct the signal using SphereFRI
[t, p, a] = SphereFRI(f, K);

% Error
[~, idx1] = sort(t);
[~, idx2] = sort(S(:, 2));

e = [(t(idx1)-S(idx2,2)) (p(idx1)-S(idx2,3)) (a(idx1)-S(idx2,1))]

sum(abs(e).^2)