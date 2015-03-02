function F = RandomSamplingMatrix(n, L)

% Generate points whose angle is uniformly distributed on 
% the unit sphere
x = rand(3, n) - 0.5;

% Get the angles (watch out, Matlab does weird things here...)n
[phi, th, ~] = cart2sph(x(1, :), x(2, :), x(3, :));
th = pi/2 - th;

i = 0;
F = zeros(n, (L+1)^2);
for l = 0:L
    for m = -l:l
        i = i + 1;
        F(:, i) = SphericalHarmonic(l, m, phi, th);       
    end
end
