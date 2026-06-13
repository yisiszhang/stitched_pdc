function compile_mex()
%COMPILE_MEX Compile optional MEX for MatrixCompletion (speed-up).
%
%   This is optional. If compilation fails, the package can still run, but
%   some matrix completion routines may be slower.
%
%   Usage:
%       setup;
%       compile_mex;

root = fileparts(mfilename('fullpath'));
mcdir = fullfile(root, '..', 'external', 'MatrixCompletion');
if exist(mcdir, 'dir') ~= 7
    error('MatrixCompletion directory not found.');
end

cwd = pwd;
cleanup = onCleanup(@() cd(cwd));
cd(mcdir);

if exist('compile.m', 'file') == 2
    compile();
else
    error('external/MatrixCompletion/compile.m not found.');
end
end
