function S = nnm_completion_complex(S, miss_mat)
%NNM_COMPLETION_COMPLEX Nuclear-norm completion for complex Hermitian spectra.
%
%   S = nnm_completion_complex(S, miss_mat) fills missing entries of a
%   complex Hermitian spectral matrix S(ω) using nuclear-norm matrix
%   completion on the real 2K-by-2K embedding.
%
%   Requires third-party `MatrixCompletion` in external/MatrixCompletion.

[~, K, nf] = size(S);

% Observed mask in complex space -> duplicated mask in real embedding
obs = ~miss_mat;
mask = [obs obs; obs obs];

lambda_tol = 10;
tol = 1e-8;
iter = 100;

for i = 1:nf
    X = S(:,:,i);
    X = [real(X) -imag(X); imag(X) real(X)];
    [Xhat, ~] = MatrixCompletion(X, mask, iter, 'nuclear', lambda_tol, tol, 0);
    S(:,:,i) = Xhat(1:K,1:K) + 1i * Xhat(K+1:2*K,1:K);
end
end
