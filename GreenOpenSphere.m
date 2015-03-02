function G = GreenOpenSphere(r, r_s, k, L)

norm_r    = sqrt(sum(r.^2));
norm_r_s  = sqrt(sum(r_s.^2));

cos_theta = r_s'*r./norm_r./norm_r_s;
cos_theta(cos_theta >  1) =  1;
cos_theta(cos_theta < -1) = -1;

G = zeros(size(cos_theta));

for l = 0:L
    P_lm = legendre(l, cos_theta);
    j_l  = SphericalBesselNaive(l, k*norm_r, 'j');
    h_l  = SphericalHankelNaive(l, k*norm_r_s);
    
    G = G + j_l .* h_l * (2*l+1) .* P_lm(1, :);
end

G = G * 1i*k/(4*pi);