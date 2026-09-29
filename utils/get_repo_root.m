function repoRoot = get_repo_root()
%GET_REPO_ROOT Top-level folder of the ETHOS_Simulations git repository.
%
%   PURPOSE:
%   Returns the git repository root so scripts can use it as
%   CONFIG.working_dir instead of a hard-coded, machine-specific path. All
%   data trees (EthosExports, RayStationFiles, SimulationResults, ...) live
%   under this folder.
%
%   INPUTS:
%   (none)
%
%   OUTPUTS:
%   repoRoot - char, absolute path of the folder that contains '.git'.
%
%   ALGORITHM:
%   1. Start from the folder this file lives in (utils/).
%   2. Walk up one parent at a time until a folder containing '.git' is found.
%      ('.git' is a folder in a normal clone and a file in a git worktree;
%      exist() finds both.) No git executable is needed.
%   3. Error if the filesystem root is reached without finding one.
%
%   EXAMPLE:
%       CONFIG.working_dir = get_repo_root();
%
%   DEPENDENCIES: none
%
%   See also: get_default_config

    repoRoot = fileparts(mfilename('fullpath'));
    while exist(fullfile(repoRoot, '.git'), 'file') == 0
        parentDir = fileparts(repoRoot);
        if strcmp(parentDir, repoRoot)
            error('get_repo_root:NotFound', ...
                'No .git folder found above %s. Is this a git checkout?', ...
                fileparts(mfilename('fullpath')));
        end
        repoRoot = parentDir;
    end
end
