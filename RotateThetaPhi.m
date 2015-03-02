function [t_rot, p_rot] = RotateThetaPhi(t, p, a, b, g, invert)

% This is a bit of a roundtrip, but at least I'm sure it works

[x, y, z] = ssht_s2c(t, p, 1);

X = [x(:)'; y(:)'; z(:)'];
R = RotationMatrix(a, b, g);

if invert
    R = inv(R);
end

X_rot = R * X;
[t_rot, p_rot, ~] = ssht_c2s(X_rot(1,:)', X_rot(2,:)', X_rot(3,:)');
