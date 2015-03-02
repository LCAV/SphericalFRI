function f_m = SelectEms(f, m, L, l_min)

if nargin == 3
    l_min = 0;
end

l_m = (max(abs(m), l_min)):L;

% Extract all spectral elements at the current m
f_m = f(l_m.^2 + l_m + m + 1);
