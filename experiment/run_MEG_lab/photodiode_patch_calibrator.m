% photodiode_patch_calibrator.m
%
% Purpose:
%   Open a standalone Psychtoolbox screen for interactively choosing the
%   photodiode patch position and size. The accepted values are printed as
%   directly reusable MATLAB assignments:
%
%     photodiodePatchPositionPx = [x y];
%     photodiodePatchSizePx = [width height];
%
% Usage example (inside experiment/):
%   run('photodiode_patch_calibrator.m');
%
% Usage example (MATLAB batch is not appropriate because this tool requires
% interactive keyboard input and a visible participant display):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment');
%   run('photodiode_patch_calibrator.m');
%
% Controls:
%   Left/Right arrows - move patch horizontally.
%   Up/Down arrows    - move patch vertically.
%   A/D               - decrease/increase patch width.
%   S/W               - decrease/increase patch height.
%   Shift             - use the fast step while held with another key.
%   Return/Enter      - accept values, print them, and exit.
%   Escape            - cancel without printing accepted assignments.
%
% Editable inputs:
%   initialPositionPx - [x y] offset from the PTB window top-left.
%   initialSizePx     - [width height] in pixels.
%   normalStepPx      - movement/resize step without Shift.
%   fastStepPx        - movement/resize step while Shift is held.
%
% Outputs:
%   On Return/Enter, the caller workspace retains:
%     photodiodePatchPositionPx
%     photodiodePatchSizePx
%   The same assignments are printed in the Command Window.
%
% Assumptions:
%   - The highest-numbered Psychtoolbox screen is the participant display.
%   - Coordinates use the same top-left-relative convention as the runtime
%     photodiode patch configuration.
%   - This is a layout utility, not a timing test; PTB sync tests are skipped
%     temporarily and their previous preference is restored on exit.

%% Editable starting values and interaction steps
% Data flow: editable defaults -> bounded patch state -> interactive updates.
clearvars;
clc;

initialPositionPx = [0 0];
initialSizePx = [50 50];
normalStepPx = 1;
fastStepPx = 10;

%% Validate utility settings before opening a display
% Reject malformed values early so the PTB window is never opened with an
% undefined interaction state.
local_validate_pair(initialPositionPx, 'initialPositionPx', 0);
local_validate_pair(initialSizePx, 'initialSizePx', 1);
local_validate_positive_integer(normalStepPx, 'normalStepPx');
local_validate_positive_integer(fastStepPx, 'fastStepPx');

initialPositionPx = round(double(initialPositionPx(:)'));
initialSizePx = round(double(initialSizePx(:)'));
normalStepPx = round(double(normalStepPx));
fastStepPx = round(double(fastStepPx));

%% Open the participant display and interactively adjust the patch
% Data flow: keyboard state -> bounded position/size -> patch rectangle ->
% display flip. The patch is deliberately drawn after the instructions.
window = [];
oldSkipSyncTests = [];
try
    try
        sca;
    catch
        % No existing PTB window is required for this standalone utility.
    end

    KbName('UnifyKeyNames');
    key = local_resolve_key_codes();

    oldSkipSyncTests = Screen('Preference', 'SkipSyncTests', 1);
    screens = Screen('Screens');
    screenNumber = max(screens);
    backgroundColor = [0 0 0];
    patchColor = [255 255 255];
    textColor = [190 190 190];
    statusColor = [255 255 0];

    [window, windowRect] = Screen('OpenWindow', screenNumber, backgroundColor);
    HideCursor;
    Screen('TextSize', window, 22);

    windowWidthPx = round(windowRect(3) - windowRect(1));
    windowHeightPx = round(windowRect(4) - windowRect(2));
    sizePx = min(initialSizePx, [windowWidthPx windowHeightPx]);
    positionPx = local_clamp_position( ...
        initialPositionPx, sizePx, windowWidthPx, windowHeightPx);

    instructions = sprintf([ ...
        'Arrow keys: move | A/D: width -/+ | S/W: height -/+\n' ...
        'Hold Shift for %d px steps | Return: accept | Escape: cancel'], ...
        fastStepPx);
    accepted = false;
    cancelled = false;

    % Discard the key that may have launched the script from the editor.
    KbReleaseWait;

    while ~(accepted || cancelled)
        Screen('FillRect', window, backgroundColor);
        DrawFormattedText(window, instructions, 'center', ...
            windowRect(4) - 105, textColor);
        statusText = sprintf('position = [%d %d] px    size = [%d %d] px', ...
            positionPx(1), positionPx(2), sizePx(1), sizePx(2));
        DrawFormattedText(window, statusText, 'center', ...
            windowRect(4) - 50, statusColor);

        % Match the runtime convention exactly: [left top right bottom],
        % with position measured from the PTB window's top-left corner.
        patchLeftTop = double(windowRect(1:2)) + double(positionPx);
        patchRect = [patchLeftTop, patchLeftTop + double(sizePx)];
        Screen('FillRect', window, patchColor, patchRect);
        Screen('Flip', window);

        [keyDown, ~, keyCode] = KbCheck;
        if ~keyDown
            WaitSecs(0.01);
            continue;
        end

        if local_any_key_pressed(keyCode, key.escape)
            cancelled = true;
            continue;
        end
        if local_any_key_pressed(keyCode, key.accept)
            accepted = true;
            continue;
        end

        shiftHeld = local_any_key_pressed(keyCode, key.shift);
        if shiftHeld
            stepPx = fastStepPx;
        else
            stepPx = normalStepPx;
        end

        if local_any_key_pressed(keyCode, key.left)
            positionPx(1) = positionPx(1) - stepPx;
        elseif local_any_key_pressed(keyCode, key.right)
            positionPx(1) = positionPx(1) + stepPx;
        elseif local_any_key_pressed(keyCode, key.up)
            positionPx(2) = positionPx(2) - stepPx;
        elseif local_any_key_pressed(keyCode, key.down)
            positionPx(2) = positionPx(2) + stepPx;
        elseif local_any_key_pressed(keyCode, key.widthDecrease)
            sizePx(1) = sizePx(1) - stepPx;
        elseif local_any_key_pressed(keyCode, key.widthIncrease)
            sizePx(1) = sizePx(1) + stepPx;
        elseif local_any_key_pressed(keyCode, key.heightDecrease)
            sizePx(2) = sizePx(2) - stepPx;
        elseif local_any_key_pressed(keyCode, key.heightIncrease)
            sizePx(2) = sizePx(2) + stepPx;
        end

        sizePx = max([1 1], min(sizePx, [windowWidthPx windowHeightPx]));
        positionPx = local_clamp_position( ...
            positionPx, sizePx, windowWidthPx, windowHeightPx);

        % A short delay provides predictable key-repeat behavior while a key
        % is held without making single-pixel adjustments feel sluggish.
        WaitSecs(0.04);
    end

    local_cleanup_calibrator(window, oldSkipSyncTests);
    window = [];

    %% Print and retain accepted values
    % Data flow: final bounded state -> runtime-compatible workspace values.
    if accepted
        photodiodePatchPositionPx = positionPx;
        photodiodePatchSizePx = sizePx;
        fprintf('\nPhotodiode patch selection accepted:\n');
        fprintf('photodiodePatchPositionPx = [%d %d];\n', ...
            photodiodePatchPositionPx(1), photodiodePatchPositionPx(2));
        fprintf('photodiodePatchSizePx = [%d %d];\n\n', ...
            photodiodePatchSizePx(1), photodiodePatchSizePx(2));
    else
        fprintf('\nPhotodiode patch selection cancelled; no values accepted.\n\n');
    end
catch ME
    local_cleanup_calibrator(window, oldSkipSyncTests);
    rethrow(ME);
end

%% Local utility functions
function local_validate_pair(value, variableName, minimumValue)
% LOCAL_VALIDATE_PAIR Validate a finite two-element integer pixel vector.
%
% Inputs:
%   value        - candidate [first second] numeric vector.
%   variableName - name used in validation errors.
%   minimumValue - inclusive minimum for both vector elements.
if ~isnumeric(value) || numel(value) ~= 2 || ...
        any(~isfinite(value)) || any(value < minimumValue) || ...
        any(value ~= round(value))
    error('%s must be a two-element integer vector with values >= %d.', ...
        variableName, minimumValue);
end
end

function local_validate_positive_integer(value, variableName)
% LOCAL_VALIDATE_POSITIVE_INTEGER Validate one positive integer pixel step.
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value < 1 || value ~= round(value)
    error('%s must be a positive integer.', variableName);
end
end

function key = local_resolve_key_codes()
% LOCAL_RESOLVE_KEY_CODES Resolve all keyboard controls after key unification.
%
% Output:
%   key - struct of scalar or vector Psychtoolbox key codes.
key = struct();
key.left = KbName('LeftArrow');
key.right = KbName('RightArrow');
key.up = KbName('UpArrow');
key.down = KbName('DownArrow');
key.widthDecrease = KbName('a');
key.widthIncrease = KbName('d');
key.heightDecrease = KbName('s');
key.heightIncrease = KbName('w');
key.escape = KbName('ESCAPE');
key.shift = [KbName('LeftShift'), KbName('RightShift')];

% Return is the main keyboard Enter key. Include keypad Enter when the
% current platform exposes a distinct key name.
key.accept = KbName('Return');
try
    keypadEnter = KbName('ENTER');
    key.accept = unique([key.accept, keypadEnter]);
catch
    % Main Return remains available when keypad Enter has no mapping.
end
end

function pressed = local_any_key_pressed(keyCode, candidateCodes)
% LOCAL_ANY_KEY_PRESSED Test safely whether any candidate code is down.
candidateCodes = candidateCodes(isfinite(candidateCodes) & candidateCodes > 0);
candidateCodes = candidateCodes(candidateCodes <= numel(keyCode));
pressed = ~isempty(candidateCodes) && any(keyCode(candidateCodes));
end

function positionPx = local_clamp_position( ...
        positionPx, sizePx, windowWidthPx, windowHeightPx)
% LOCAL_CLAMP_POSITION Keep the complete patch inside the PTB window.
%
% Inputs:
%   positionPx                - requested top-left-relative [x y].
%   sizePx                    - current [width height].
%   windowWidthPx/HeightPx    - display extent in pixels.
% Output:
%   positionPx                - bounded integer [x y].
maxPositionPx = max([0 0], ...
    [windowWidthPx windowHeightPx] - double(sizePx));
positionPx = min(max(round(double(positionPx(:)')), [0 0]), maxPositionPx);
end

function local_cleanup_calibrator(window, oldSkipSyncTests)
% LOCAL_CLEANUP_CALIBRATOR Restore PTB state after normal exit or failure.
%
% Inputs:
%   window           - current PTB window handle, or [] before window setup.
%   oldSkipSyncTests - previous SkipSyncTests preference, or [] if unknown.
if ~isempty(window)
    try
        Screen('Close', window);
    catch
        % Continue restoring remaining global PTB state.
    end
end
try
    Priority(0);
catch
end
try
    ShowCursor;
catch
end
if ~isempty(oldSkipSyncTests)
    try
        Screen('Preference', 'SkipSyncTests', oldSkipSyncTests);
    catch
        % Preference restoration is best-effort if PTB is unavailable.
    end
end
end
