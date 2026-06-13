function [Theta, W, info] = graphicalLassoDiag(S, rho, maxIt, tol)
%GRAPHICALLASSODIAG Original graphical lasso update with iteration diagnostics.

p = size(S,1);
if nargin < 4, tol = 1e-6; end
if nargin < 3, maxIt = 1e2; end

W = S + rho * eye(p);
W_old = W;
i = 0;
converged = false;
maxShootingIterations = 0;
nShootingMaxed = 0;

while i < maxIt
    i = i + 1;
    for j = p:-1:1
        jminus = setdiff(1:p,j);
        [V,D] = eig(W(jminus,jminus));
        d = diag(D);
        X = V * diag(sqrt(d)) * V';
        Y = V * diag(1./sqrt(d)) * V' * S(jminus,j);
        [b, shootInfo] = local_lasso_shooting(X, Y, rho, maxIt, tol);
        maxShootingIterations = max(maxShootingIterations, shootInfo.iterations);
        nShootingMaxed = nShootingMaxed + double(shootInfo.hit_maxiter);
        W(jminus,j) = W(jminus,jminus) * b;
        W(j,jminus) = W(jminus,j)';
    end
    if norm(W - W_old, 1) < tol
        converged = true;
        break;
    end
    W_old = W;
end

Theta = W^-1;
info.iterations = i;
info.max_iterations = maxIt;
info.hit_maxiter = (i == maxIt) && ~converged;
info.converged = converged;
info.max_shooting_iterations = maxShootingIterations;
info.n_shooting_maxed = nShootingMaxed;
end


function [b, info] = local_lasso_shooting(X, Y, lambda, maxIt, tol)
if nargin < 4, tol = 1e-6; end
if nargin < 3, maxIt = 1e2; end

[n,p] = size(X);
if p > n
    b = zeros(p,1);
else
    b = X \ Y;
end
b_old = b;
i = 0;
converged = false;

XTX = X' * X;
XTY = X' * Y;

while i < maxIt
    i = i + 1;
    for j = 1:p
        jminus = setdiff(1:p,j);
        S0 = XTX(j,jminus) * b(jminus) - XTY(j);
        if S0 > lambda
            b(j) = (lambda - S0) / norm(X(:,j),2)^2;
        elseif S0 < -lambda
            b(j) = -(lambda + S0) / norm(X(:,j),2)^2;
        else
            b(j) = 0;
        end
    end
    if norm(b - b_old, 1) < tol
        converged = true;
        break;
    end
    b_old = b;
end

info.iterations = i;
info.hit_maxiter = (i == maxIt) && ~converged;
end
