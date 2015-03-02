function CRLB = CramerRaoBound(theta, phi, sigma2, t0, p0, a0, L, posonly)

DfDt0 = zeros(size(theta(:)));
DfDp0 = zeros(size(theta(:)));
DfDa0 = zeros(size(theta(:)));

for l = 0:(L-1)
    for m = -l:l
        Y_lm_n  = SphericalHarmonic(l, m, theta, phi);
        Y_lm1_0 = SphericalHarmonic(l, m+1, t0, p0);
        Y_lm_0  = SphericalHarmonic(l, m, t0, p0);
        
        DfDa0 = DfDa0 + conj(Y_lm_0) * Y_lm_n;
        DfDt0 = DfDt0 + a0 * (m*cot(t0)*conj(Y_lm_0) ...
              + sqrt((l-m)*(l+m+1)) * exp(1i*p0) ...
             .* conj(Y_lm1_0)) * Y_lm_n;
        DfDp0 = DfDp0 + a0 * (-1i*m) * conj(Y_lm_0) * Y_lm_n;
    end
end

if posonly
    Df = [DfDt0 DfDp0].';
else
    Df = [DfDt0 DfDp0 DfDa0].';
end


% Df   = [DfDt0 DfDp0].';
I    = 1/sigma2 * (Df * Df');
CRLB = inv(I);
