function cm = PolypartCoeffs(l, m, L)
% Note that this doesn't contain any normalization. Normalization should be
% appended typically in the creation of spherical harmonics. On the other
% hand, it includes (-1)^m (Condon-Shortley style).

c0 = zeros(L+1, 1);
c0(end-l:end) = LegendreCoeffs(l);

% Derivative matrix
D  = diag(L:-1:1, -1);

% Coefficients
cm = D^abs(m) * c0;

% To satisfy Condon-Shortley phase convention
cm = (-1)^m * cm;

    
