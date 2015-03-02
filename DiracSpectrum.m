function f = DiracSpectrum(S, L)

i = 0;
f = zeros(L^2, 1);
for l = 0:L-1
    for m = -l:l
        i = i + 1;
        Y    = conj(SphericalHarmonic(l, m, S(:, 2), S(:, 3)));
        f(i) = sum(S(:, 1) .* Y);
    end
end
