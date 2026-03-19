% Detailed batch version (line-by-line comments)

% Absolute path to the one-dot script you want to run.
scriptPath = '/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/one-dot/PE_simulation_diff_1Dot.m';

% Participants to process.
participants = [73 72 71];

% The three new PCR strategies to test.
strategies = {'baseline_pcr_border','ridge_full_autocorr', 'ridge_tapered_autocorr', 'ar1_prewhite_ridge'};

% Loop participants.
for iP = 1:numel(participants)
    % Loop strategies.
    for iS = 1:numel(strategies)
        % Run one participant + one strategy with explicit parameters.
        run_one_dot_with_explicit_params(participants(iP), strategies{iS}, scriptPath);
    end
end

function run_one_dot_with_explicit_params(pid, strategy, scriptPath)
    % Subject ID to load (MovDot_SubXX.mat).
    participantNumber = pid;

    % Select PCR branch of the pipeline.
    dRSAtypeToRun = 'PCR';

    % Run only nondeviant condition (change to {'nondeviant','deviant'} if needed).
    inputConditions = {'nondeviant','deviant'};

    % 2 = force full rerun and overwrite matching existing outputs.
    existingResultsAction = 2;

    % 0 = keep normal console progress output.
    suppressDispText = 0;

    % Sampling rate used by the simulation script.
    plotSampleRateHz = 120;

    % Container for parameter overrides.
    params = struct();

    % PCR-specific parameter container.
    params.PCR = struct();

    % Choose one of: ridge_full_autocorr / ridge_tapered_autocorr / ar1_prewhite_ridge.
    params.PCR.RegressStrategy = strategy;

    % Ridge strength: tuned default for AR1 strategy vs other ridge strategies.
    if strcmp(strategy, 'ar1_prewhite_ridge')
        params.PCR.RidgeLambdaFactor = 0.04;   % recommended for AR1+ridge
    else
        params.PCR.RidgeLambdaFactor = 0.06;   % recommended for full/tapered ridge
    end

    % Width of lag-distance Gaussian (used by tapered penalty profile).
    params.PCR.TaperSigmaFactor = 0.33;

    % How much extra penalty to add to far lags in tapered mode (>=0).
    params.PCR.AutoPenaltyStrength = 2;

    % 1 = normalize mean autocorr penalty per row to keep edge/center comparable.
    params.PCR.NormalizeAutoPenaltyPerRow = 1;

    % 1 = z-score predictors before ridge solve.
    params.PCR.StandardizePredictors = 1;

    % Stability clip for pooled AR(1) coefficient in ar1_prewhite_ridge.
    params.PCR.AR1Clip = 0.99;

    % Execute the pipeline.
    run(scriptPath);
end
