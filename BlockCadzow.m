function Z = BlockCadzow(Z, K, L, reality)

n_iter = 100;

for iter = 1:n_iter
    
    ind = 0;
    for m = -(L-1-K):(L-1-K)

        % First average along the diagonals

        Y = Z(ind+1:ind+(L-abs(m)-K), :);
        [r, n] = size(Y);
        for d = 0:(r+n-2)
            M = MaskAdiag(r, n, d);
            Y(logical(M)) = mean(Y(logical(M)));
        end
        Z(ind+1:ind+(L-abs(m)-K), :) = Y;
        ind = ind + (L-abs(m)-K);
    end
    
    % Then do SVD thresholding

    [U, S, V] = svd(Z);
    S(:, K+1) = 0;
    Z = U*S*V';
end    
