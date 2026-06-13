function [Sred, fred] = ibl_reduce_frequency_grid(S, f, targetN)
%IBL_REDUCE_FREQUENCY_GRID Reduce frequency dimension by contiguous-bin averaging.
%
% The reduced frequency axis is returned on an exactly uniform grid so that
% downstream Wilson factorization, which checks for even spacing, remains valid.

nf = numel(f);
if nargin < 3 || isempty(targetN) || targetN <= 0 || nf <= targetN
    Sred = S;
    fred = f(:);
    return;
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

Sred = zeros(sz(1), sz(2), targetN, 'like', S);
for k = 1:targetN
    idx = edges(k):(edges(k+1)-1);
    Sred(:,:,k) = mean(S(:,:,idx), 3, 'omitnan');
end

% Wilson factorization expects an evenly spaced frequency axis. Use a
% uniform grid spanning the reduced band instead of block-center averages,
% which become slightly irregular when bins have unequal widths.
fred = linspace(f(1), f(end), targetN).';

if wasVector
    Sred = squeeze(Sred);
end
end
