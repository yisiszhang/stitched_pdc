function run_selftest()
%RUN_SELFTEST Minimal smoke test to catch missing dependencies/paths.

fprintf('Running SNPDC self-test...\n');
setup();

% Check key functions
req = {'reconstruct_inversepsd', 'stitched_pdc', 'nonparam_pdc_H', 'arsim', 'sfactorization_wilson'};
for i = 1:numel(req)
    assert(exist(req{i}, 'file') == 2, 'Missing function on path: %s', req{i});
end

% Quick VAR(1) tiny run (no figures)
A = [0.5 -0.1 0.1; -0.4 0.5 0; -0.1 0.2 0.3];
C = eye(3)*0.1;
w = zeros(3,1);
n = 2000; ndisc = 200;

recset = { [1 2], [2 3] };
x = cell(1,2);
for u = 1:2
    v = arsim(w, A, C, n, ndisc);
    x{u} = v(:, recset{u});
end

params.fs = 1;
params.win = bartlett(16);
params.nov = 10;
params.nfft = 256;
params.lambda = 0;
params.method = 'nnm';

[~, ~, Sreg, f] = reconstruct_inversepsd(x, recset, params);
PDC = nonparam_pdc_H(Sreg, f);
assert(all(isfinite(PDC(:))), 'PDC contains non-finite values.');

fprintf('OK.\n');
end
