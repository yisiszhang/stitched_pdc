function S = maxdet_completion_complex(S, miss_mat)
%MAXDET_COMPLETION_COMPLEX Max-determinant completion for complex Hermitian spectra.
%
%   S = maxdet_completion_complex(S, miss_mat) imputes missing entries
%   of S(ω) by maximizing log det of the completed covariance in the real
%   2K-by-2K embedding, subject to matching observed entries.
%
%   Note: This routine is intended for small/medium K due to the use of
%   generic optimization (fminunc).

[~,n,nf] = size(S);

ind = [[miss_mat miss_mat]; [miss_mat miss_mat]];

options = optimset('GradObj', 'on', 'MaxIter', 1000); 

for i = 1:nf
    x0 = S(:,:,i);
    x0 = [[real(x0) -imag(x0)]; [imag(x0) real(x0)]];
    [x, cost] = fminunc(@(x)(costFunction(x, ind)), x0(:), options);    
    x = reshape(x, n*2, n*2);
    S(:,:,i) = x(1:n, 1:n) + sqrt(-1) * x(n+1:n*2, 1:n);
end


function [cost, grad] = costFunction(x, ind)
    l = length(x);
    x = reshape(x, sqrt(l), sqrt(l));
    cost = - det(x);

    grad = -adjoint(x)';
    grad(ind==0) = 0;
    grad = grad(:);
end

end
