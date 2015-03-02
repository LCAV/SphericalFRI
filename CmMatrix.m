function C_m = CmMatrix(m, L_max, l_min)

if nargin == 2
    l_min = 0;
end

l_min = max(abs(m), l_min);
l_m = l_min:L_max;
C_m = zeros(length(l_m));

for l = l_m
    c_lm = PolypartCoeffs(l, abs(m), L_max);
    N_lm = sqrt((2*l+1)*factorial(l-abs(m))/factorial(l+abs(m))/4/pi);

    if (m < 0)
        c_lm = (-1)^m * c_lm;
    end

    C_m(l-l_min+1, :) = N_lm * c_lm((abs(m)+1):(abs(m)+length(l_m)));
end
