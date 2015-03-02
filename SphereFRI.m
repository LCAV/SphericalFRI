function [theta, phi, alpha] = SphereFRI(f, K, doBlockCadzow, doRelaxation)

if nargin < 3
    doBlockCadzow = false;
    doRelaxation = false;
elseif nargin < 4
    doRelaxation = false;
end

% Spectrum length
L = round(sqrt(length(f)));

% Get the parts of X*A*U that we can get, and construct the annihilating
% matrix
Z = zeros(L-K + sum(1:L-K-1)*2, K+1);
r = 0;
for m = -(L-1):(L-1)
    C_m = CmMatrix(m, L-1);
    f_m = SelectEms(f, m, L-1); % Extract spectral elements at current m
    
    % Invert to get the linear combination of cosine powers in X*A*u_m
    y = C_m \ f_m;
    
    % Construct the 'annihilation' matrix rows (Toeplitz by parts).
    % Non-zero number of iterations will occur for m < L-K. Note that in
    % the paper we have a Hankel-by-parts formulation.
    for i = 1:(L-abs(m)-K)
        r = r + 1;
        Z(r, :) = y(i:i+K);
    end
end

if doBlockCadzow
    Z = BlockCadzow(Z, K, L, 0);
end

% Find the nullspace of Z (use SVD for robustness, annihilation filter is
% the last column of V)
[~, ~, V] = svd(Z);

% Find the colatitudes
ctheta = real(roots(V(:, end))); % Matlab roots assumes that the first 
                                 % entry corresponds to the largest power
theta  = acos(ctheta);

% Now construct the matrix X, and invert to find U(:, m==1). Note that
% below for Au_1, the largest power is (L-1), since for m==1, there is no 
% power L.
X    = bsxfun(@power, repmat(ctheta', L, 1), ((L-1):-1:0)');
Au_0 = (CmMatrix(0, L-1)*X) \ SelectEms(f, 0, L-1);
Au_1 = (CmMatrix(1, L-1)*X(2:end, :)) \ SelectEms(f, 1, L-1);

phi = mod(-angle(Au_1 ./ Au_0), 2*pi);

% Find the peak amplitudes
alpha = Au_0;

if doRelaxation
    [theta, phi, alpha] = RelaxDiracs(theta, phi, alpha, f);
end

