function h = SphericalHankelNaive(n, z)

h = sqrt(pi/2./z) .* besselh(n + 0.5, z);
