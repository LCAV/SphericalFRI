function S = RandomDiracs(K, cpx)

S = rand(K, 3); % (mags, thetas, phis)
S = S * diag([1, pi, 2*pi]);

if cpx
    S(:, 1) = S(:, 1) + rand(K, 1) * 1i;
end

