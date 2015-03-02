function Y = SphericalHarmonic(l, m, theta, phi)

theta = theta(:);
phi   = phi(:);

if abs(m) > l
    Y = zeros(size(theta));
    return;
end

% Use abs(m) to avoid division by 0

c = PolypartCoeffs(l, m, l);
X = bsxfun(@power, repmat(cos(theta'), l+1, 1), (l:-1:0)');
N = sqrt((2*l+1)*factorial(l-abs(m))/factorial(l+abs(m))/4/pi);

Y = N * (c' * X)' .* sin(theta).^abs(m) .* exp(1i*m*phi);

if (m < 0)
    Y = (-1)^m * Y; % no need for conjugation because the argument of exp up 
                    % there has 'm', not 'abs(m)'
end
