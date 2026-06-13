function setup(varargin)
%SETUP  Setup paths and optional MEX compilation for snpdc
%
% Usage:
%   setup                % add paths only
%   setup compile        % add paths and compile optional MEX
%
% Notes:
%   Some optional acceleration (MatrixCompletion) requires compiling
%   a MEX file locally due to macOS Gatekeeper restrictions.


fprintf('[snpdc] Setting up paths...\n');

root = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(root,'src')));

% External dependencies
extDir = fullfile(root,'external');
if exist(extDir,'dir')
    addpath(genpath(extDir));
end

% Note: external/sfactorization_wilson.m is from FieldTrip (GPL) and is added via genpath above.


% ------------------------------------------------------------
% Optional MEX compilation
% ------------------------------------------------------------
doCompile = (nargin > 0) && strcmpi(varargin{1}, 'compile');

mexName = 'MinimizeSquaredDiffbudget';
mexFile = [mexName, '.', mexext];
mexPath = which(mexName);

if isempty(mexPath)
    fprintf('[snpdc] Optional MEX "%s" not found.\n', mexFile);

    if doCompile
        fprintf('[snpdc] Attempting to compile MEX...\n');
        try
            mcDir = fullfile(root,'external','MatrixCompletion');
            if ~exist(mcDir,'dir')
                error('MatrixCompletion directory not found.');
            end
            cur = pwd;
            cd(mcDir);
            srcFile = [mexName, '.cpp'];
            if exist(srcFile, 'file') ~= 2
                srcFile = [mexName, '.c'];
            end
            mex(srcFile);
            cd(cur);

            fprintf('[snpdc] MEX compilation successful.\n');
        catch ME
            fprintf('[snpdc] MEX compilation failed.\n');
            fprintf('[snpdc] Reason: %s\n', ME.message);
            fprintf('[snpdc] Continuing without MEX (MATLAB-only mode).\n');
        end
    else
        fprintf('[snpdc] Matrix completion will run in MATLAB-only mode (slower).\n');
        fprintf('[snpdc] To enable acceleration, run:\n');
        fprintf('         setup compile\n');
    end
else
    fprintf('[snpdc] Found MEX: %s\n', mexPath);
end

fprintf('[snpdc] Setup complete.\n');

end
