function v = SphericalBesselNaive(n, z, kind)

if kind == 'j'
    v = sqrt(pi/2./z) .* besselj(n + 1/2, z);
    v(z == 0) = (n == 0);
elseif kind == 'y'
    v = sqrt(pi/2./z) .* bessely(n + 1/2, z);
end    


