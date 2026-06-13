function [Sout, fout, idx] = reduce_frequency_grid(S, f, targetN)
%REDUCE_FREQUENCY_GRID Downsample spectra to a manageable frequency grid.
%
%   [Sout, fout, idx] = reduce_frequency_grid(S, f, targetN)
%
% Uses contiguous-bin averaging and returns an exactly uniform frequency
% axis. The uniform axis matters because Wilson factorization checks that
% frequencies are evenly spaced.

nf = numel(f);
if nargin < 3 || isempty(targetN) || targetN <= 0 || nf <= targetN
    Sout = S;
    fout = f(:);
    idx = (1:nf).';
    return
end

edges = round(linspace(1, nf + 1, targetN + 1));
edges(end) = nf + 1;

sz = size(S);
if isvector(S)
    S = reshape(S, 1, 1, []);
    wasVector = true;
else
    wasVector = false;
end

Sout = zeros(sz(1), sz(2), targetN, 'like', S);
idx = cell(targetN, 1);
for k = 1:targetN
    idx{k} = edges(k):(edges(k+1)-1);
    Sout(:,:,k) = mean(S(:,:,idx{k}), 3, 'omitnan');
end

fout = linspace(f(1), f(end), targetN).';

if wasVector
    Sout = squeeze(Sout);
end
end
