function [x, y, z] = SphCart(t, p, r)

t = pi/2 - t;
[x, y, z] = sph2cart(p, t, r);
