function C = meacount_mat(recset, K)
%MEACOUNT_MAT Count how many times each pair of variables is co-observed.
%
%   C(i,j) is the number of blocks u such that i and j are both in recset{u}.
%
%   If K is omitted, it is inferred from the maximum index appearing in recset.

if nargin < 2 || isempty(K)
    K = 0;
    for u = 1:numel(recset)
        if ~isempty(recset{u})
            K = max(K, max(recset{u}(:)));
        end
    end
end

C = zeros(K, K);
for u = 1:numel(recset)
    idx = recset{u};
    if isempty(idx), continue; end
    idx = unique(idx(:))';
    idx = idx(idx >= 1 & idx <= K);
    if isempty(idx), continue; end
    C(idx, idx) = C(idx, idx) + 1;
end
end
