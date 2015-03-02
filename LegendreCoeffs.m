function p = LegendreCoeffs(l)

% Document somewhere all these tricks that you're using. Perhaps just by
% making a Matlab notebook.

p = zeros(l+1, 1);

for k = 0:floor(l/2)
    p(l-2*k+1) = (-1)^k * nchoosek(l, k) * nchoosek(2*l-2*k, l);
end

p = flipud(p) / 2^l;

