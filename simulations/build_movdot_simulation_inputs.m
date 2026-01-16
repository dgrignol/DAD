%% Build dot-path inputs for dRSA simulations (MoveDot1 stimuli)
% Script: build_movdot_simulation_inputs.m
%
% Purpose:
%   Load a MovDot_SubXX.mat stimulus file, split trials into deviant vs.
%   non-deviant conditions, and save per-condition 3D arrays for each dot
%   (trials × features × time) in visual degrees. Output filenames include
%   the dot colors defined by MoveDot1_experiment_vX.m.
%
% Example usage (from repo root in MATLAB):
%   addpath('simulations');
%   build_movdot_simulation_inputs('experiment/input_files/MovDot_Sub98.mat');
%
% Example usage (custom output folder and MoveDot script):
%   build_movdot_simulation_inputs('experiment/input_files/MovDot_Sub98.mat', ...
%       'OutputDir', 'simulations/input', ...
%       'MoveDotScript', 'experiment/MoveDot1_experiment_vX.m');
%
% Inputs:
%   - stimuliPath (char/string): Full or relative path to MovDot_SubXX.mat.
%       Relative paths are checked against the current working directory
%       first, then the repo root (parent of simulations).
%   - Name/value pairs:
%       'OutputDir'     : output folder for .mat files (default: simulations/input).
%                         Relative paths are anchored to the repo root unless
%                         they already exist relative to the working directory.
%       'MoveDotScript' : path to MoveDot1_experiment_vX.m for dot colors.
%                         Relative paths are resolved like stimuliPath.
%
% Outputs:
%   - Saves two .mat files (deviant, non-deviant) into OutputDir, each with:
%       dot1Paths (trials × 2 × time) in visual degrees
%       dot2Paths (trials × 2 × time) in visual degrees
%       meta (struct with source info and color labels)
%
% Key assumptions:
%   - xySeqs(...).xy columns are [x1 y1 x2 y2] in visual degrees.
%   - Deviant trials have condition > 0; non-deviant have condition == 0.
%   - All trials within a condition share the same number of frames.
function outputFiles = build_movdot_simulation_inputs(stimuliPath, varargin)
    %% Parse inputs and establish default paths
    parser = inputParser;
    parser.addRequired('stimuliPath', @(x) ischar(x) || isstring(x));

    scriptDir = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(scriptDir);
    defaultOutputDir = fullfile(scriptDir, 'input');
    defaultMoveDotScript = fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m');

    parser.addParameter('OutputDir', defaultOutputDir, @(x) ischar(x) || isstring(x));
    parser.addParameter('MoveDotScript', defaultMoveDotScript, @(x) ischar(x) || isstring(x));
    parser.parse(stimuliPath, varargin{:});

    stimuliPath = char(parser.Results.stimuliPath);
    outputDir = char(parser.Results.OutputDir);
    moveDotScript = char(parser.Results.MoveDotScript);

    %% Resolve input/output paths (repo-root aware)
    % Data flow: user-provided paths -> resolved absolute paths -> file IO.
    stimuliPath = resolve_existing_path(stimuliPath, repoRoot, scriptDir, 'file', true);
    outputDir = resolve_output_path(outputDir, repoRoot);
    moveDotScript = resolve_existing_path(moveDotScript, repoRoot, scriptDir, 'file', false);

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    %% Load stimulus paths and derive condition split
    data = load(stimuliPath);
    if ~isfield(data, 'xySeqs')
        error('Stimulus file does not contain xySeqs: %s', stimuliPath);
    end
    xyAll = data.xySeqs(:);
    conditions = [xyAll.condition];
    deviantMask = conditions > 0;
    nonDeviantMask = conditions == 0;

    %% Resolve dot colors from MoveDot1_experiment_vX.m (for filenames)
    [dot1Rgb, dot2Rgb] = read_dot_colors(moveDotScript);
    dot1Name = rgb_to_name(dot1Rgb);
    dot2Name = rgb_to_name(dot2Rgb);

    %% Build per-condition arrays and save results
    [~, baseName, ~] = fileparts(stimuliPath);
    outputFiles = struct('condition', {}, 'file', {});

    outputFiles(end+1) = save_condition( ...
        baseName, 'nondeviant', nonDeviantMask, xyAll, outputDir, ...
        dot1Rgb, dot2Rgb, dot1Name, dot2Name, stimuliPath);
    outputFiles(end+1) = save_condition( ...
        baseName, 'deviant', deviantMask, xyAll, outputDir, ...
        dot1Rgb, dot2Rgb, dot1Name, dot2Name, stimuliPath);
end

%% Helper: resolve existing file/directory paths with repo-root fallback
function resolvedPath = resolve_existing_path(rawPath, repoRoot, scriptDir, pathKind, mustExist)
    if is_absolute_path(rawPath)
        resolvedPath = rawPath;
        if mustExist && ~path_exists(resolvedPath, pathKind)
            error('Path not found: %s', rawPath);
        end
        return;
    end

    % Check current working directory first, then repo root, then script dir.
    candidates = {rawPath, fullfile(repoRoot, rawPath), fullfile(scriptDir, rawPath)};
    resolvedPath = '';
    for i = 1:numel(candidates)
        if path_exists(candidates{i}, pathKind)
            resolvedPath = candidates{i};
            break;
        end
    end

    if isempty(resolvedPath)
        if mustExist
            error('Path not found: %s (checked: %s)', rawPath, strjoin(candidates, ', '));
        end
        resolvedPath = fullfile(repoRoot, rawPath);
    end
end

%% Helper: resolve output directory path (repo-root default for relative)
function outputDir = resolve_output_path(outputDir, repoRoot)
    if is_absolute_path(outputDir)
        return;
    end
    if exist(outputDir, 'dir')
        return;
    end

    % For relative directories, prefer the repo root to avoid cwd surprises.
    outputDir = fullfile(repoRoot, outputDir);
end

%% Helper: detect absolute paths (cross-platform)
function tf = is_absolute_path(pathIn)
    if ispc
        tf = ~isempty(regexp(pathIn, '^[A-Za-z]:[\\/]', 'once')) || startsWith(pathIn, '\\');
    else
        tf = startsWith(pathIn, filesep);
    end
end

%% Helper: check file/dir existence based on type
function tf = path_exists(pathIn, pathKind)
    switch lower(pathKind)
        case 'file'
            tf = isfile(pathIn);
        case 'dir'
            tf = exist(pathIn, 'dir') ~= 0;
        otherwise
            error('Unknown path kind: %s', pathKind);
    end
end

%% Helper: save a condition-specific file
function entry = save_condition(baseName, conditionLabel, mask, xyAll, outputDir, ...
        dot1Rgb, dot2Rgb, dot1Name, dot2Name, stimuliPath)
    trials = xyAll(mask);
    if isempty(trials)
        warning('No trials found for condition: %s', conditionLabel);
        entry.condition = conditionLabel;
        entry.file = '';
        return;
    end

    % Validate frame consistency to keep a clean trials × features × time cube.
    nFrames = arrayfun(@(s) size(s.xy, 1), trials);
    uniqueFrames = unique(nFrames);
    if numel(uniqueFrames) ~= 1
        error('Condition %s has variable frame counts: %s', ...
            conditionLabel, mat2str(uniqueFrames));
    end

    nTrials = numel(trials);
    nTime = uniqueFrames(1);
    dot1GreenPaths = zeros(nTrials, 2, nTime);
    dot2YellowPaths = zeros(nTrials, 2, nTime);

    % Populate dot1/dot2 arrays: each feature slice is x or y.
    for iTrial = 1:nTrials
        xy = trials(iTrial).xy;
        dot1GreenPaths(iTrial, 1, :) = reshape(xy(:, 1), 1, 1, []);
        dot1GreenPaths(iTrial, 2, :) = reshape(xy(:, 2), 1, 1, []);
        dot2YellowPaths(iTrial, 1, :) = reshape(xy(:, 3), 1, 1, []);
        dot2YellowPaths(iTrial, 2, :) = reshape(xy(:, 4), 1, 1, []);
    end

    meta = struct();
    meta.sourceFile = stimuliPath;
    meta.condition = conditionLabel;
    meta.units = 'visual_degrees';
    meta.features = {'x', 'y'};
    meta.dot1Rgb = dot1Rgb;
    meta.dot2Rgb = dot2Rgb;
    meta.dot1Name = dot1Name;
    meta.dot2Name = dot2Name;
    meta.deviantDefinition = 'condition > 0 (non-deviant = condition == 0)';

    outName = sprintf('%s_%s.mat', baseName, conditionLabel);
    outPath = fullfile(outputDir, outName);
    save(outPath, 'dot1GreenPaths', 'dot2YellowPaths', 'meta');

    entry.condition = conditionLabel;
    entry.file = outPath;
end

%% Helper: read Conf.DotColor1/2 from MoveDot1_experiment_vX.m
function [dot1Rgb, dot2Rgb] = read_dot_colors(moveDotScript)
    % Defaults if the script is missing or parsing fails.
    dot1Rgb = [0 255 0];
    dot2Rgb = [242 223 0];

    if ~isfile(moveDotScript)
        warning('MoveDot script not found; using default dot colors.');
        return;
    end

    lines = read_text_lines(moveDotScript);
    dot1Rgb = parse_rgb_assignment(lines, 'Conf.DotColor1', dot1Rgb);
    dot2Rgb = parse_rgb_assignment(lines, 'Conf.DotColor2', dot2Rgb);
end

%% Helper: parse RGB assignment from a line list
function rgb = parse_rgb_assignment(lines, token, fallback)
    rgb = fallback;
    expr = token + "\\s*=\\s*\\[([^\\]]+)\\]";
    for i = 1:numel(lines)
        line = strtrim(lines(i));
        matches = regexp(line, expr, 'tokens', 'once');
        if ~isempty(matches)
            values = sscanf(matches{1}, '%f');
            if numel(values) >= 3
                rgb = values(1:3)';
                return;
            end
        end
    end
end

%% Helper: read text lines (compatibility with older MATLAB)
function lines = read_text_lines(filePath)
    % Use fileread/splitlines to avoid dependency on readlines.
    rawText = fileread(filePath);
    lines = splitlines(string(rawText));
end

%% Helper: map RGB (0-255) to a simple color name
function name = rgb_to_name(rgb)
    colorTable = struct( ...
        'name', {'green', 'yellow', 'red', 'blue', 'white', 'black', 'gray', 'orange'}, ...
        'rgb', {[0 255 0], [255 255 0], [255 0 0], [0 0 255], [255 255 255], [0 0 0], [128 128 128], [255 165 0]});

    dists = zeros(1, numel(colorTable));
    for i = 1:numel(colorTable)
        diff = double(rgb) - double(colorTable(i).rgb);
        dists(i) = sqrt(sum(diff.^2));
    end
    [~, idx] = min(dists);
    name = colorTable(idx).name;
end
