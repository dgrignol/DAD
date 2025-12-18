# Decreasing dot speed by 25% (to 3/4) and lengthening trials

Edit these parameters, then regenerate inputs.

- `experiment/lib/Config.m`
  - `dotSpeedDegPerSec`: multiply by `0.75` (e.g., `4.973 -> 3.73`). `dotSpeedDegPerFrame` updates automatically.
  - `trialDuration`: multiply by `4/3` (e.g., `2 -> 2.67`) so trials last longer at the slower speed.
  - Per-stimulus path durations (scale each by `4/3`):
    - If using `Config.likelihood`: change `likelihood.pathDuration` (e.g., `2.0 -> 2.67`).
    - If using `Config.path_duration`: scale every value in `path_duration.pathDuration`.
    - If using `Config.path_duration_norm`: scale every value in `path_duration_norm.pathDuration`.
  - Leave `frameFrequency` unchanged.

- Optional (only if you want occlusion displacement to stay similar despite slower dots):
  - `experiment/CreateInputFiles_v10.m`: scale `Catch.OcclDuration` by `4/3` (e.g., `0.2 -> 0.27`) so `Catch.OcclDistance` stays comparable when dot speed drops.

Regenerate files after edits (per subject):
1) Run `experiment/stimuli_generation_v5.m` to recreate `input_files/MovDot_SubXX.mat` with the slower speed/longer trials.
2) Run `experiment/CreateInputFiles_v10.m` to rebuild `input_files/SubXX_TrialStruct.mat` (catch timing uses the new `fps`/`dpf`/frame counts).
