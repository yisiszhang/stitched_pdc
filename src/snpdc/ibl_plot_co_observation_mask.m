function out = ibl_plot_co_observation_mask(qc, orderMode)
%IBL_PLOT_CO_OBSERVATION_MASK Plot post-QC co-observation support across sessions.
%
%   out = ibl_plot_co_observation_mask(qc)
%   out = ibl_plot_co_observation_mask(qc, orderMode)
%
% Input
%   qc   struct returned by ibl_build_pca_qc_graph
%   orderMode  'original', 'cluster', or 'symamd' (default: 'cluster')
%
% Output
%   out  struct containing support matrices used for plotting

if nargin < 2 || isempty(orderMode)
    orderMode = 'cluster';
end

areaNames = string(qc.area_names(:));
nAreas = numel(areaNames);
nSessions = numel(qc.qualifying_sessions);

pairMask = false(nAreas, nAreas);
pairCount = zeros(nAreas, nAreas);

for k = 1:nSessions
    sessionAreas = string(qc.qualifying_sessions(k).area_names(:));
    [tf, idx] = ismember(sessionAreas, areaNames);
    idx = idx(tf);
    if numel(idx) < 2
        continue;
    end
    pairMask(idx, idx) = true;
    pairCount(idx, idx) = pairCount(idx, idx) + 1;
end

pairFrac = pairCount ./ max(nSessions, 1);
[perm, resolvedMode] = local_reorder(pairCount, orderMode);

areaNamesPlot = areaNames(perm);
pairMaskPlot = pairMask(perm, perm);
pairFracPlot = pairFrac(perm, perm);
pairCountPlot = pairCount(perm, perm);

figure('Name', 'IBL Co-Observation Mask', 'Color', 'w');
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(pairMaskPlot);
axis image;
colormap(gca, gray(2));
colorbar('Ticks', [0 1], 'TickLabels', {'no', 'yes'});
title(sprintf('Co-observation mask (%d sessions, %s order)', nSessions, resolvedMode));
set(gca, 'XTick', 1:nAreas, 'XTickLabel', areaNamesPlot, ...
    'YTick', 1:nAreas, 'YTickLabel', areaNamesPlot, 'XTickLabelRotation', 90);

nexttile;
imagesc(pairFracPlot, [0 1]);
axis image;
colorbar;
title(sprintf('Co-observation fraction (%s order)', resolvedMode));
set(gca, 'XTick', 1:nAreas, 'XTickLabel', areaNamesPlot, ...
    'YTick', 1:nAreas, 'YTickLabel', areaNamesPlot, 'XTickLabelRotation', 90);

sgtitle(sprintf('Qualified-session support for %d post-QC areas', nAreas));

out.area_names = areaNames;
out.area_names_plot = areaNamesPlot;
out.pair_mask = pairMask;
out.pair_count = pairCount;
out.pair_fraction = pairFrac;
out.pair_mask_plot = pairMaskPlot;
out.pair_count_plot = pairCountPlot;
out.pair_fraction_plot = pairFracPlot;
out.n_sessions = nSessions;
out.permutation = perm;
out.order_mode = resolvedMode;
end


function [perm, resolvedMode] = local_reorder(pairCount, orderMode)
nAreas = size(pairCount, 1);
perm = (1:nAreas)';
resolvedMode = lower(string(orderMode));

switch resolvedMode
    case "original"
        return;
    case "symamd"
        W = pairCount + pairCount.';
        perm = symamd(sparse(W));
        perm = perm(:);
        return;
    case "cluster"
        if exist('linkage', 'file') == 2 && nAreas >= 3
            rowVecs = pairCount;
            rowNorm = sqrt(sum(rowVecs.^2, 2));
            rowNorm(rowNorm == 0) = 1;
            rowVecs = rowVecs ./ rowNorm;
            sim = rowVecs * rowVecs.';
            sim = max(min(sim, 1), -1);
            dist = 1 - sim;
            dist = (dist + dist.') / 2;
            dist(1:nAreas+1:end) = 0;
            y = local_squareform_upper(dist);
            Z = linkage(y, 'average');
            if exist('optimalleaforder', 'file') == 2
                try
                    leafOrd = optimalleaforder(Z, y);
                    perm = leafOrd(:);
                    return;
                catch
                    % Fall back below if optimalleaforder is unavailable or fails.
                end
            end
            resolvedMode = "symamd";
            W = pairCount + pairCount.';
            perm = symamd(sparse(W));
            perm = perm(:);
            return;
        end
        resolvedMode = "symamd";
        W = pairCount + pairCount.';
        perm = symamd(sparse(W));
        perm = perm(:);
        return;
    otherwise
        error('Unknown orderMode: %s', char(orderMode));
end
end


function y = local_squareform_upper(D)
n = size(D,1);
y = zeros(n*(n-1)/2, 1);
t = 1;
for i = 1:n-1
    for j = i+1:n
        y(t) = D(i,j);
        t = t + 1;
    end
end
end
