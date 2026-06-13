function v = orient_loading(v)
%ORIENT_LOADING Fix arbitrary complex PC phase/sign convention.
%
% The largest-magnitude loading is rotated to be real positive.

if isempty(v) || ~any(isfinite(v))
    return
end
[~, idx] = max(abs(v));
anchor = v(idx);
if anchor == 0 || ~isfinite(anchor)
    return
end
v = v .* exp(-1i * angle(anchor));
if real(v(idx)) < 0
    v = -v;
end
if max(abs(imag(v))) < 1e-12
    v = real(v);
end
end
