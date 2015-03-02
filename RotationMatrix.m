function R = RotationMatrix(alpha, beta, gamma)

mtx = @(x) [cos(x) -sin(x); sin(x) cos(x)];

Rz1 = eye(3);
Ry  = eye(3);
Rz2 = eye(3);

Rz1(1:2, 1:2)    = mtx(alpha);
Ry([1,3], [1,3]) = mtx( beta);
Rz2(1:2, 1:2)    = mtx(gamma);

R = Rz1 * Ry * Rz2;