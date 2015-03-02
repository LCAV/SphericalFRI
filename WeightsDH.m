function a = WeightsDH(L)

j = 0:(2*L-1);
j = j + 1/2;
theta = pi*j / 2 / (L);
ct    = cos(theta);

P = zeros(2*L);
for l = 0:(2*L-1)
    P_lm = legendre(l, ct);
    P(l+1, :) = P_lm(1, :);
end

b = zeros(2*L, 1);
b(1) = sqrt(2);
a = P \ b;