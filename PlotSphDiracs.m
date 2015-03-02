function PlotSphDiracs(S, scale, d, nPoints)

[x, y, z] = sphere(nPoints);
[p, t, r] = cart2sph(x, y, z);

for i = 1:size(S, 1)
    r((t-S(i,2)).^2 + abs(p-S(i,3)).^2 < d^2) = 1+S(i,1)*scale;
end

[x, y, z] = sph2cart(p, t, r);

surf(x, y, z, ones(size(x)));

shading interp;
material dull;
lighting phong;
light;
axis off equal;
