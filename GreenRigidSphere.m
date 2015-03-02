function G = GreenRigidSphere(r, r_s, k, L)

norm_r    = sqrt(sum(r.^2));
norm_r_s  = sqrt(sum(r_s.^2));

cos_theta = r_s'*r./norm_r./norm_r_s;
cos_theta(cos_theta >  1) =  1;
cos_theta(cos_theta < -1) = -1;

G = zeros(size(cos_theta));

for l = 0:(L-1)
    P_lm    = legendre(l, cos_theta);
    h_prime = 1/2 * (SphericalHankelNaive(l-1, k*norm_r) ...
                     - 1./(k*norm_r).*(SphericalHankelNaive(l, k*norm_r) ...
                     + (k*norm_r).*SphericalHankelNaive(l+1, k*norm_r)));
    
    b_l  = 1i ./ h_prime ./ (k*norm_r).^2;

    G = G + b_l .* SphericalHankelNaive(l, k*norm_r_s) * (2*l+1) .* P_lm(1, :);
end
G = G * 1i*k/(4*pi);