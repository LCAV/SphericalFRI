function [theta, phi, alpha] = SphereFRI_SN(f, K, L)

% Here, once you estimate the theta's, you could round them to the nearest
% grid point and construct the matrix X according to this knowledge. That
% would improve the performance.

% Compute the bandwidth
L_prime = round(sqrt(length(f)));

% Get the parts of X*A*U that we can get, and construct the annihilating
% matrix
Z = zeros((L_prime - L - K)*(L_prime - L - K + 1), K+1);
r = 0;
for m = [(-(L_prime-1):(-L)) (L:(L_prime-1))]
    C_m = CmMatrix(m, L_prime-1);
    f_m = SelectEms(f, m, L_prime-1); % Extract spectral elements at current m
    
    % Invert to get the linear combination of cosine powers in X*A*u_m
    y = C_m \ f_m;
    
    % Construct the 'annihilation' matrix rows (Toeplitz by parts).
    % Non-zero number of iterations will occur for m <= L-K (when l_min=0).
    for i = 1:(L_prime-abs(m)-K)
        r = r + 1;
        Z(r, :) = y(i:i+K);
    end
end

% Find the nullspace of Z (use SVD for robustness, annihilation filter is
% the last column of V)
[~, ~, V] = svd(Z);

% Find the colatitudes
ctheta = real(roots(V(:, end)));
theta  = acos(ctheta);

% Now construct the matrix X, and invert to find U(:, m==1). Note that
% below the largest power is (L-1), since for m==1, there is no power L.
% t = unique(t);
% dtt = distance(t', theta');
% [~, idxidx] = min(dtt);
% theta = reshape(t(idxidx), numel(theta), 1);

X    = bsxfun(@power, repmat(ctheta', L_prime, 1), ((L_prime-1):-1:0)');
Au_0 = (CmMatrix(L, L_prime-1)*X(L+1:end, :)) \ SelectEms(f, L, L_prime-1);
Au_1 = (CmMatrix(L+1, L_prime-1)*X(L+2:end, :)) \ SelectEms(f, L+1, L_prime-1);

phi = mod(-angle(Au_1 ./ Au_0), 2*pi);

% Find the peak amplitudes (why does this work in l_min=0 case? Probably
% because K is still smaller than whatever L we're left with).
u_1   = sin(theta).^(L+1) .* exp(-1i*1*phi*(L+1));
alpha = Au_1 ./ u_1;

