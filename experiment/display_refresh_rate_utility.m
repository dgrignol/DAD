function report = display_refresh_rate_utility( ...
        requestedRefreshHz, requestedWidthPx, requestedHeightPx, ...
        keepRequestedModeOnPass, screenNumber)
% DISPLAY_REFRESH_RATE_UTILITY Configure and validate a display refresh rate.
%
% Purpose:
%   Independently test the display path used by the experiment without
%   launching, modifying, or depending on the experiment runtime. The utility:
%     1) measures the current Psychtoolbox flip cadence;
%     2) requests the desired resolution and refresh rate from the OS/GPU;
%     3) measures the resulting flip cadence over many real flips;
%     4) queries the DATAPixx video input status when available;
%     5) prints a compact PASS/FAIL summary.
%
% Usage examples (inside experiment/):
%   display_refresh_rate_utility
%   report = display_refresh_rate_utility(120);
%   report = display_refresh_rate_utility(120, 1920, 1080, true);
%
% To test the requested mode but restore the original mode afterwards:
%   report = display_refresh_rate_utility(120, 1920, 1080, false);
%
% To test a specific Psychtoolbox screen instead of the highest-numbered one:
%   report = display_refresh_rate_utility(120, 1920, 1080, true, 1);
%
% Inputs:
%   requestedRefreshHz        - desired GPU output refresh in Hz; default 120.
%   requestedWidthPx          - desired display width; default 1920.
%   requestedHeightPx         - desired display height; default 1080.
%   keepRequestedModeOnPass   - true keeps the requested mode after PASS;
%                               false restores the original mode; default true.
%   screenNumber              - Psychtoolbox screen id; default max(Screen).
%
% Output:
%   report - struct containing the original/requested modes, pre-request and
%            post-request PTB measurements, DATAPixx status, PASS/FAIL checks,
%            and whether the requested mode was left active.
%
% Example Command Window summary:
%   Pre-request:      100.000 Hz
%   Requested:        120.000 Hz
%   PTB measured:     119.985 Hz
%   Median IFI:         8.334 ms
%   DATAPixx input:   120.000 Hz
%   Result:           PASS
%
% Key assumptions and safety behavior:
%   - Native PROPixx 1920x1080 at 120 Hz is driven by the GPU video mode; this
%     utility does not enable a PROPixx high-speed sequencer.
%   - The requested mode must be exposed by Screen('Resolutions'), which may
%     require prior VPixx EDID and operating-system display configuration.
%   - Psychtoolbox synchronization tests remain enabled throughout timing
%     measurements. An existing PTB window is treated as an error.
%   - A failed request, failed validation, or unexpected error restores the
%     original display mode. A successful request is kept only when requested.
%   - DATAPixx verification is required for PASS because the intended MEG
%     acquisition path includes the DATAPixx controller.

%% Resolve defaults and validate caller inputs
% Data flow: optional arguments -> validated target mode and display id.
if nargin < 1 || isempty(requestedRefreshHz)
    requestedRefreshHz = 120;
end
if nargin < 2 || isempty(requestedWidthPx)
    requestedWidthPx = 1920;
end
if nargin < 3 || isempty(requestedHeightPx)
    requestedHeightPx = 1080;
end
if nargin < 4 || isempty(keepRequestedModeOnPass)
    keepRequestedModeOnPass = true;
end

local_validate_positive_scalar(requestedRefreshHz, 'requestedRefreshHz');
local_validate_positive_integer(requestedWidthPx, 'requestedWidthPx');
local_validate_positive_integer(requestedHeightPx, 'requestedHeightPx');
if ~(islogical(keepRequestedModeOnPass) || ...
        (isnumeric(keepRequestedModeOnPass) && ...
        isscalar(keepRequestedModeOnPass) && ...
        isfinite(keepRequestedModeOnPass) && ...
        ismember(double(keepRequestedModeOnPass), [0 1])))
    error('keepRequestedModeOnPass must be logical or numeric 0/1.');
end

requestedRefreshHz = double(requestedRefreshHz);
requestedWidthPx = round(double(requestedWidthPx));
requestedHeightPx = round(double(requestedHeightPx));
keepRequestedModeOnPass = logical(keepRequestedModeOnPass);

AssertOpenGL;
screens = Screen('Screens');
if isempty(screens)
    error('Psychtoolbox did not report any display screens.');
end
if nargin < 5 || isempty(screenNumber)
    screenNumber = max(screens);
else
    local_validate_nonnegative_integer(screenNumber, 'screenNumber');
    screenNumber = round(double(screenNumber));
    if ~ismember(screenNumber, screens)
        error('screenNumber %d is not in Screen(''Screens'') = %s.', ...
            screenNumber, mat2str(screens));
    end
end

% This utility changes a global display mode and must never disrupt a running
% experiment or another PTB onscreen window.
openWindows = Screen('Windows');
if ~isempty(openWindows)
    error(['Close all Psychtoolbox windows before running this utility. ' ...
        'Open window handles: %s'], mat2str(openWindows));
end

%% Configure validation thresholds and report state
% Data flow: requested rate -> acceptance thresholds -> final PASS/FAIL checks.
measurementFlipCount = 360;
warmupFlipCount = 30;
refreshToleranceHz = 0.50;
intervalToleranceMs = 0.75;
minimumIntervalsWithinTolerancePct = 99.0;

report = struct();
report.screenNumber = screenNumber;
report.requestedMode = struct( ...
    'width', requestedWidthPx, ...
    'height', requestedHeightPx, ...
    'hz', requestedRefreshHz);
report.keepRequestedModeOnPass = keepRequestedModeOnPass;
report.originalMode = struct();
report.availableModes = struct([]);
report.requestedModeAvailable = false;
report.preRequest = local_empty_timing_stats();
report.postRequest = local_empty_timing_stats();
report.postRequestOsMode = struct();
report.dataPixx = local_empty_datapixx_status();
report.checks = struct();
report.result = 'ERROR';
report.requestedModeKept = false;
report.originalModeRestored = false;

oldSkipSyncTests = [];
originalModeKnown = false;
modeRequestAttempted = false;
modeCurrentlyRequested = false;

try
    %% Snapshot the original mode and measure the unmodified display
    % Data flow: active OS/GPU mode -> PTB flip timestamps -> pre-request Hz.
    oldSkipSyncTests = Screen('Preference', 'SkipSyncTests', 0);
    report.originalMode = Screen('Resolution', screenNumber);
    originalModeKnown = true;
    report.availableModes = Screen('Resolutions', screenNumber);
    report.requestedModeAvailable = local_mode_is_available( ...
        report.availableModes, requestedWidthPx, requestedHeightPx, ...
        requestedRefreshHz, refreshToleranceHz);

    fprintf('\nDisplay refresh utility: screen %d\n', screenNumber);
    fprintf('Original OS mode: %s\n', local_format_mode(report.originalMode));
    fprintf('Measuring pre-request refresh (%d flips)...\n', ...
        measurementFlipCount);
    report.preRequest = local_measure_flip_timing( ...
        screenNumber, measurementFlipCount, warmupFlipCount, ...
        requestedRefreshHz, intervalToleranceMs);
    fprintf('Pre-request:     %9.3f Hz (median IFI %.3f ms)\n', ...
        report.preRequest.measuredRefreshHz, ...
        report.preRequest.medianIfiMs);

    %% Verify availability and request the target GPU video mode
    % Data flow: Screen('Resolutions') + target -> OS/GPU mode request.
    fprintf('Requested:       %9.3f Hz (%dx%d)\n', ...
        requestedRefreshHz, requestedWidthPx, requestedHeightPx);
    if ~report.requestedModeAvailable
        availableRates = local_rates_for_resolution( ...
            report.availableModes, requestedWidthPx, requestedHeightPx);
        if isempty(availableRates)
            availableText = 'none';
        else
            availableText = mat2str(availableRates, 6);
        end
        fprintf(['Requested mode is not exposed by Screen(''Resolutions'').\n' ...
            'Available rates for %dx%d: %s Hz\n'], ...
            requestedWidthPx, requestedHeightPx, availableText);
        report.result = 'FAIL';
        report.checks = local_build_checks( ...
            report, requestedRefreshHz, requestedWidthPx, ...
            requestedHeightPx, refreshToleranceHz, ...
            minimumIntervalsWithinTolerancePct);
        local_print_summary(report);
        local_restore_skip_sync_tests(oldSkipSyncTests);
        return;
    end

    modeRequestAttempted = true;
    local_set_resolution(screenNumber, requestedWidthPx, ...
        requestedHeightPx, requestedRefreshHz, report.originalMode);
    modeCurrentlyRequested = true;

    % Allow the OS, GPU, DATAPixx, and projector video link to relock before
    % opening the validation window.
    WaitSecs(2.0);
    report.postRequestOsMode = Screen('Resolution', screenNumber);

    %% Measure the requested mode and query the DATAPixx input
    % Data flow: requested OS/GPU mode -> PTB VBL timestamps + controller
    % video status -> independent timing and signal-path validation.
    fprintf('Measuring post-request refresh (%d flips)...\n', ...
        measurementFlipCount);
    report.postRequest = local_measure_flip_timing( ...
        screenNumber, measurementFlipCount, warmupFlipCount, ...
        requestedRefreshHz, intervalToleranceMs);
    report.dataPixx = local_query_datapixx_status();

    %% Evaluate all rate, resolution, stability, and controller checks
    % A PASS requires agreement between the requested mode, PTB-measured
    % cadence, and the DATAPixx input seen on the acquisition video path.
    report.checks = local_build_checks( ...
        report, requestedRefreshHz, requestedWidthPx, ...
        requestedHeightPx, refreshToleranceHz, ...
        minimumIntervalsWithinTolerancePct);
    checkValues = struct2cell(report.checks);
    if all(cellfun(@(value) isequal(value, true), checkValues))
        report.result = 'PASS';
    else
        report.result = 'FAIL';
    end

    local_print_summary(report);

    %% Keep a successful requested mode or restore the original mode
    % Data flow: validation result + caller policy -> final desktop mode.
    if strcmp(report.result, 'PASS') && keepRequestedModeOnPass
        report.requestedModeKept = true;
        fprintf('Display mode:    requested mode kept active.\n\n');
    else
        local_restore_original_mode(screenNumber, report.originalMode);
        modeCurrentlyRequested = false;
        report.originalModeRestored = true;
        if strcmp(report.result, 'PASS')
            fprintf('Display mode:    original mode restored by request.\n\n');
        else
            fprintf('Display mode:    original mode restored after FAIL.\n\n');
        end
    end

    local_restore_skip_sync_tests(oldSkipSyncTests);
catch ME
    % Close any partially opened PTB state before attempting a resolution
    % restoration. The local measurement function normally closes its own
    % window, but this is a final safeguard for unexpected PTB failures.
    try
        Screen('CloseAll');
    catch
    end
    try
        Priority(0);
    catch
    end
    try
        ShowCursor;
    catch
    end

    if originalModeKnown && (modeRequestAttempted || modeCurrentlyRequested)
        try
            local_restore_original_mode(screenNumber, report.originalMode);
        catch restoreError
            restoreMessage = sprintf( ...
                'Could not restore the original display mode: %s', ...
                restoreError.message);
            warning('displayRefreshUtility:ModeRestoreFailed', '%s', ...
                restoreMessage);
        end
    end
    local_restore_skip_sync_tests(oldSkipSyncTests);
    rethrow(ME);
end
end

%% Local measurement and hardware-query functions
function stats = local_measure_flip_timing( ...
        screenNumber, measurementFlipCount, warmupFlipCount, ...
        targetRefreshHz, intervalToleranceMs)
% LOCAL_MEASURE_FLIP_TIMING Measure VBL cadence over scheduled PTB flips.
%
% Inputs:
%   screenNumber          - target Psychtoolbox screen id.
%   measurementFlipCount  - number of retained inter-flip intervals.
%   warmupFlipCount       - scheduled flips discarded before measurement.
%   targetRefreshHz       - rate used only for interval-stability scoring.
%   intervalToleranceMs   - absolute per-interval target tolerance.
% Output:
%   stats - measured refresh, IFI distribution, and missed-frame diagnostics.
window = [];
stats = local_empty_timing_stats();
try
    [window, windowRect] = Screen('OpenWindow', screenNumber, [0 0 0]);
    HideCursor;
    Priority(MaxPriority(window));

    measuredIfiSec = double(Screen('GetFlipInterval', window));
    if ~isfinite(measuredIfiSec) || measuredIfiSec <= 0
        error('Screen(''GetFlipInterval'') returned invalid IFI %.9f s.', ...
            measuredIfiSec);
    end

    % Alternate a very dark full-screen value so every retained flip carries
    % new image content without producing a bright flickering test display.
    vbl = Screen('Flip', window);
    nRetainedTimestamps = measurementFlipCount + 1;
    retainedVbl = nan(nRetainedTimestamps, 1);
    retainedMissed = nan(nRetainedTimestamps, 1);
    retainedIndex = 0;
    totalScheduledFlips = warmupFlipCount + nRetainedTimestamps;

    for iFlip = 1:totalScheduledFlips
        darkLevel = mod(iFlip, 2);
        Screen('FillRect', window, [darkLevel darkLevel darkLevel], ...
            windowRect);
        requestedWhen = vbl + 0.5 * measuredIfiSec;
        [vbl, ~, ~, missed] = Screen('Flip', window, requestedWhen);

        if iFlip > warmupFlipCount
            retainedIndex = retainedIndex + 1;
            retainedVbl(retainedIndex) = vbl;
            retainedMissed(retainedIndex) = missed;
        end
    end

    intervalMs = diff(retainedVbl) * 1000;
    intervalMs = intervalMs(isfinite(intervalMs) & intervalMs > 0);
    if numel(intervalMs) ~= measurementFlipCount
        error('Expected %d valid IFIs but measured %d.', ...
            measurementFlipCount, numel(intervalMs));
    end

    sortedIntervals = sort(intervalMs);
    stats.getFlipIntervalMs = measuredIfiSec * 1000;
    stats.medianIfiMs = median(intervalMs);
    stats.meanIfiMs = mean(intervalMs);
    stats.stdIfiMs = std(intervalMs);
    stats.p01IfiMs = local_percentile_sorted(sortedIntervals, 0.01);
    stats.p99IfiMs = local_percentile_sorted(sortedIntervals, 0.99);
    stats.measuredRefreshHz = 1000 / stats.medianIfiMs;
    stats.nIntervals = numel(intervalMs);
    stats.positiveMissedCount = sum( ...
        isfinite(retainedMissed) & retainedMissed > 0);

    targetIfiMs = 1000 / targetRefreshHz;
    stats.intervalsWithinTargetTolerancePct = 100 * mean( ...
        abs(intervalMs - targetIfiMs) <= intervalToleranceMs);
    stats.longIntervalCount = sum(intervalMs > 1.5 * targetIfiMs);

    Screen('Close', window);
    window = [];
    Priority(0);
    ShowCursor;
catch ME
    if ~isempty(window)
        try
            Screen('Close', window);
        catch
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
    rethrow(ME);
end
end

function statusInfo = local_query_datapixx_status()
% LOCAL_QUERY_DATAPIXX_STATUS Read the controller's current video input mode.
%
% Output:
%   statusInfo - availability, receiving state, input resolution/frequency,
%                and any diagnostic message if the query cannot be completed.
statusInfo = local_empty_datapixx_status();
dataPixxOpened = false;
if exist('Datapixx', 'file') == 0
    statusInfo.message = 'Datapixx toolbox is not on the MATLAB path.';
    return;
end

try
    Datapixx('Open');
    dataPixxOpened = true;
    Datapixx('RegWrRd');
    rawStatus = Datapixx('GetVideoStatus');

    statusInfo.available = true;
    statusInfo.receivingVideo = logical( ...
        local_struct_numeric_field(rawStatus, 'receivingVideo', 0));
    statusInfo.widthPx = local_struct_numeric_field( ...
        rawStatus, 'horizontalResolution', NaN);
    statusInfo.heightPx = local_struct_numeric_field( ...
        rawStatus, 'verticalResolution', NaN);
    statusInfo.verticalFrequencyHz = local_struct_numeric_field( ...
        rawStatus, 'verticalFrequency', NaN);
    statusInfo.message = '';
catch ME
    statusInfo.message = ME.message;
end

if dataPixxOpened
    try
        Datapixx('Close');
    catch closeError
        if isempty(statusInfo.message)
            statusInfo.message = sprintf( ...
                'Video status read succeeded, but close failed: %s', ...
                closeError.message);
        end
    end
end
end

%% Local validation, mode handling, and reporting functions
function checks = local_build_checks( ...
        report, targetHz, targetWidth, targetHeight, ...
        refreshToleranceHz, minimumIntervalsWithinTolerancePct)
% LOCAL_BUILD_CHECKS Convert independent measurements into PASS criteria.
%
% A requested mode is valid only if the OS exposes it, the post-request mode
% has the expected size, PTB measures the expected cadence and stability, and
% DATAPixx reports the same incoming video mode.
checks = struct();
checks.requestedModeAvailable = logical(report.requestedModeAvailable);

osMode = report.postRequestOsMode;
checks.osResolutionMatches = ...
    local_mode_dimension_matches(osMode, targetWidth, targetHeight);
checks.osRefreshMatches = local_mode_refresh_matches( ...
    osMode, targetHz, refreshToleranceHz);

ptb = report.postRequest;
checks.ptbRefreshMatches = isfinite(ptb.measuredRefreshHz) && ...
    abs(ptb.measuredRefreshHz - targetHz) <= refreshToleranceHz;
checks.ptbIntervalsStable = ...
    isfinite(ptb.intervalsWithinTargetTolerancePct) && ...
    ptb.intervalsWithinTargetTolerancePct >= ...
    minimumIntervalsWithinTolerancePct;
checks.ptbNoLongIntervals = isfinite(ptb.longIntervalCount) && ...
    ptb.longIntervalCount == 0;

dp = report.dataPixx;
checks.dataPixxAvailable = logical(dp.available);
checks.dataPixxReceivingVideo = logical(dp.receivingVideo);
checks.dataPixxResolutionMatches = ...
    isfinite(dp.widthPx) && isfinite(dp.heightPx) && ...
    round(dp.widthPx) == targetWidth && round(dp.heightPx) == targetHeight;
checks.dataPixxRefreshMatches = isfinite(dp.verticalFrequencyHz) && ...
    abs(dp.verticalFrequencyHz - targetHz) <= refreshToleranceHz;
end

function local_print_summary(report)
% LOCAL_PRINT_SUMMARY Print the requested compact pre/post comparison.
fprintf('\nRefresh validation summary\n');
fprintf('Pre-request:     %9s\n', ...
    local_format_hz(report.preRequest.measuredRefreshHz));
fprintf('Requested:       %9.3f Hz\n', report.requestedMode.hz);
fprintf('PTB measured:    %9s\n', ...
    local_format_hz(report.postRequest.measuredRefreshHz));
fprintf('Median IFI:      %9s\n', ...
    local_format_ms(report.postRequest.medianIfiMs));

if report.dataPixx.available
    fprintf('DATAPixx input:  %9s (%gx%g, receiving=%d)\n', ...
        local_format_hz(report.dataPixx.verticalFrequencyHz), ...
        report.dataPixx.widthPx, report.dataPixx.heightPx, ...
        double(report.dataPixx.receivingVideo));
else
    fprintf('DATAPixx input:  unavailable\n');
    if ~isempty(report.dataPixx.message)
        fprintf('DATAPixx note:   %s\n', report.dataPixx.message);
    end
end
fprintf('Result:          %s\n', report.result);

if strcmp(report.result, 'FAIL') && ~isempty(fieldnames(report.checks))
    failedNames = fieldnames(report.checks);
    failedMask = false(size(failedNames));
    for iField = 1:numel(failedNames)
        failedMask(iField) = ~isequal( ...
            report.checks.(failedNames{iField}), true);
    end
    failedNames = failedNames(failedMask);
    if ~isempty(failedNames)
        fprintf('Failed checks:   %s\n', strjoin(failedNames, ', '));
    end
end
end

function available = local_mode_is_available( ...
        modes, targetWidth, targetHeight, targetHz, toleranceHz)
% LOCAL_MODE_IS_AVAILABLE Test whether Screen exposes the requested mode.
available = false;
for iMode = 1:numel(modes)
    mode = modes(iMode);
    if local_mode_dimension_matches(mode, targetWidth, targetHeight) && ...
            local_mode_refresh_matches(mode, targetHz, toleranceHz)
        available = true;
        return;
    end
end
end

function matches = local_mode_dimension_matches(mode, targetWidth, targetHeight)
% LOCAL_MODE_DIMENSION_MATCHES Compare a Screen mode struct to target size.
matches = isstruct(mode) && isfield(mode, 'width') && ...
    isfield(mode, 'height') && isfinite(double(mode.width)) && ...
    isfinite(double(mode.height)) && ...
    round(double(mode.width)) == targetWidth && ...
    round(double(mode.height)) == targetHeight;
end

function matches = local_mode_refresh_matches(mode, targetHz, toleranceHz)
% LOCAL_MODE_REFRESH_MATCHES Compare a reported mode rate with tolerance.
matches = isstruct(mode) && isfield(mode, 'hz') && ...
    isfinite(double(mode.hz)) && double(mode.hz) > 0 && ...
    abs(double(mode.hz) - targetHz) <= toleranceHz;
end

function rates = local_rates_for_resolution(modes, targetWidth, targetHeight)
% LOCAL_RATES_FOR_RESOLUTION List refresh rates exposed for one resolution.
rates = [];
for iMode = 1:numel(modes)
    mode = modes(iMode);
    if local_mode_dimension_matches(mode, targetWidth, targetHeight) && ...
            isfield(mode, 'hz') && isfinite(double(mode.hz))
        rates(end + 1) = double(mode.hz); %#ok<AGROW>
    end
end
rates = unique(rates);
end

function local_set_resolution( ...
        screenNumber, widthPx, heightPx, refreshHz, originalMode)
% LOCAL_SET_RESOLUTION Request a supported OS/GPU display mode.
%
% Preserve the original color depth when Screen reports one, because reducing
% pixel depth can disable alpha blending and alter experiment rendering.
if isfield(originalMode, 'pixelSize') && ...
        isfinite(double(originalMode.pixelSize)) && ...
        double(originalMode.pixelSize) > 0
    Screen('Resolution', screenNumber, widthPx, heightPx, refreshHz, ...
        double(originalMode.pixelSize));
else
    Screen('Resolution', screenNumber, widthPx, heightPx, refreshHz);
end
end

function local_restore_original_mode(screenNumber, originalMode)
% LOCAL_RESTORE_ORIGINAL_MODE Restore the complete pre-test display mode.
if ~isstruct(originalMode) || ...
        ~all(isfield(originalMode, {'width', 'height', 'hz'}))
    error('Original display mode is incomplete and cannot be restored.');
end

widthPx = double(originalMode.width);
heightPx = double(originalMode.height);
refreshHz = double(originalMode.hz);
if isfield(originalMode, 'pixelSize') && ...
        isfinite(double(originalMode.pixelSize)) && ...
        double(originalMode.pixelSize) > 0
    Screen('Resolution', screenNumber, widthPx, heightPx, refreshHz, ...
        double(originalMode.pixelSize));
else
    Screen('Resolution', screenNumber, widthPx, heightPx, refreshHz);
end
WaitSecs(1.0);
end

function textValue = local_format_mode(mode)
% LOCAL_FORMAT_MODE Render a Screen resolution struct for human inspection.
if ~isstruct(mode) || ...
        ~all(isfield(mode, {'width', 'height', 'hz'}))
    textValue = 'unavailable';
    return;
end
textValue = sprintf('%gx%g at %.3f Hz', ...
    double(mode.width), double(mode.height), double(mode.hz));
end

function textValue = local_format_hz(value)
% LOCAL_FORMAT_HZ Render a finite refresh value or an unavailable marker.
if isfinite(value)
    textValue = sprintf('%.3f Hz', value);
else
    textValue = 'unavailable';
end
end

function textValue = local_format_ms(value)
% LOCAL_FORMAT_MS Render a finite millisecond value or unavailable marker.
if isfinite(value)
    textValue = sprintf('%.3f ms', value);
else
    textValue = 'unavailable';
end
end

function value = local_struct_numeric_field(rawStruct, fieldName, fallback)
% LOCAL_STRUCT_NUMERIC_FIELD Read one finite scalar hardware-status value.
value = fallback;
if isstruct(rawStruct) && isfield(rawStruct, fieldName)
    candidate = double(rawStruct.(fieldName));
    if isscalar(candidate) && isfinite(candidate)
        value = candidate;
    end
end
end

function percentileValue = local_percentile_sorted(sortedValues, fraction)
% LOCAL_PERCENTILE_SORTED Select a bounded nearest-rank percentile.
index = max(1, min(numel(sortedValues), ...
    ceil(double(fraction) * numel(sortedValues))));
percentileValue = sortedValues(index);
end

function stats = local_empty_timing_stats()
% LOCAL_EMPTY_TIMING_STATS Create a stable schema before measurements exist.
stats = struct( ...
    'getFlipIntervalMs', NaN, ...
    'medianIfiMs', NaN, ...
    'meanIfiMs', NaN, ...
    'stdIfiMs', NaN, ...
    'p01IfiMs', NaN, ...
    'p99IfiMs', NaN, ...
    'measuredRefreshHz', NaN, ...
    'nIntervals', 0, ...
    'positiveMissedCount', NaN, ...
    'intervalsWithinTargetTolerancePct', NaN, ...
    'longIntervalCount', NaN);
end

function statusInfo = local_empty_datapixx_status()
% LOCAL_EMPTY_DATAPIXX_STATUS Create a stable controller-status schema.
statusInfo = struct( ...
    'available', false, ...
    'receivingVideo', false, ...
    'widthPx', NaN, ...
    'heightPx', NaN, ...
    'verticalFrequencyHz', NaN, ...
    'message', '');
end

function local_restore_skip_sync_tests(oldValue)
% LOCAL_RESTORE_SKIP_SYNC_TESTS Restore the caller's global PTB preference.
if ~isempty(oldValue)
    try
        Screen('Preference', 'SkipSyncTests', oldValue);
    catch restoreError
        restoreMessage = sprintf( ...
            'Could not restore SkipSyncTests preference: %s', ...
            restoreError.message);
        warning('displayRefreshUtility:PreferenceRestoreFailed', '%s', ...
            restoreMessage);
    end
end
end

function local_validate_positive_scalar(value, variableName)
% LOCAL_VALIDATE_POSITIVE_SCALAR Validate one finite value greater than zero.
if ~isnumeric(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    error('%s must be a finite numeric scalar > 0.', variableName);
end
end

function local_validate_positive_integer(value, variableName)
% LOCAL_VALIDATE_POSITIVE_INTEGER Validate one finite positive integer.
if ~isnumeric(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 1 || value ~= round(value)
    error('%s must be a finite positive integer.', variableName);
end
end

function local_validate_nonnegative_integer(value, variableName)
% LOCAL_VALIDATE_NONNEGATIVE_INTEGER Validate one finite integer >= 0.
if ~isnumeric(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 0 || value ~= round(value)
    error('%s must be a finite integer >= 0.', variableName);
end
end
