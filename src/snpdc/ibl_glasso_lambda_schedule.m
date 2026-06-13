function lambda = ibl_glasso_lambda_schedule(S, lambda0, power)
%IBL_GLASSO_LAMBDA_SCHEDULE Build a frequency-specific glasso penalty.
%
%   lambda = ibl_glasso_lambda_schedule(S, lambda0)
%   lambda = ibl_glasso_lambda_schedule(S, lambda0, power)
%
% Uses lambda(f) = lambda0 * median(abs(diag(S(:,:,f))))^power.

if nargin < 2 || isempty(lambda0)
    lambda0 = 1e-6;
end
if nargin < 3 || isempty(power)
    power = 3;
end

nf = size(S, 3);
lambda = zeros(nf, 1);
for f = 1:nf
    diagPower = abs(diag(S(:,:,f)));
    diagPower = diagPower(isfinite(diagPower));
    if isempty(diagPower)
        scale = 0;
    else
        scale = median(diagPower);
    end
    lambda(f) = lambda0 * scale.^power;
end
lambda(~isfinite(lambda) | lambda < 0) = 0;
end
