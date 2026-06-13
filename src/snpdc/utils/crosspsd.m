function [S, f] = crosspsd(X, window, noverlap, nfft, fs)
%CROSSPSD Estimate cross power spectral density matrix using Welch's method.
%
%   [S,f] = crosspsd(X, window, noverlap, nfft, fs) returns a 3-D array S of
%   size (N x N x F), where N is the number of channels (columns of X) and F
%   is the number of frequency bins returned by CPSD. S(i,j,:) is the CPSD
%   estimate between channels i and j.
%
%   This is a convenience wrapper around MATLAB's cpsd() for multichannel
%   data. Requires the Signal Processing Toolbox.

[nT, n] = size(X);
if nT < 2
    error('crosspsd:NotEnoughData', 'X must have at least 2 time points.');
end

% Precompute frequency vector using first pair.
[~, f] = cpsd(X(:,1), X(:,1), window, noverlap, nfft, fs);
nf = numel(f);

S = zeros(n, n, nf);
for i = 1:n
    for j = i:n
        Pij = cpsd(X(:,i), X(:,j), window, noverlap, nfft, fs);
        S(i,j,:) = Pij;
        if i ~= j
            S(j,i,:) = conj(Pij);
        end
    end
end
end
