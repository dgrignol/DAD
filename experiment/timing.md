# Timing Assessment for `MoveDot1_experiment_occlusion_v14_PathBandOccluder.m`

Date: 2026-04-02

## Scope and method
This assessment is based on static code inspection of:
- `experiment/MoveDot1_experiment_occlusion_v14_PathBandOccluder.m`
- `experiment/stimuli_generation_V24_PathBandOccluder.m`
- `experiment/lib/Config_stimuli_generation_V24_PathBandOccluder.m`
- `experiment/lib/Config_schedule_CreateInputV17_MoveDotV14_PathBandOccluder.m`

No live PTB timing capture was run in this step (no photodiode/MEG hardware timing trace collected here).

## Executive summary (MEG-oriented)
- The occluder centerline geometry is **precomputed offline** in stimuli generation and stored in input fields (`pathband_pre_xy`, `pathband_post_xy`, timing frames). This is good.
- In runtime, occluders are **still reconstructed and rasterized each frame** (straight-band polygon rebuilt from the polyline at draw time). This adds CPU/GPU work but is usually manageable at 120 Hz with debug overlays off.
- Event triggers are emitted in the frame loop at the configured frame indices, but the TTL write is done **before** `Screen('Flip')`, so visual-vs-trigger alignment is not guaranteed to be exactly zero-lag without external validation.
- The script currently forces `Screen('Preference','SkipSyncTests',1)`, which is not ideal for MEG precision validation.

Overall: likely usable for MEG if run on a stable 120 Hz setup with debug overlays off, but there are clear timing-hardening improvements needed before calling it timing-robust.

## Where occluder computations happen

### Precomputed (offline, not per frame)
From `stimuli_generation_V24_PathBandOccluder.m`:
- Builds pathband center polylines and occlusion metadata (`occlusion_*_frame`).
- Stores in each trial struct:
  - `pathband_pre_xy`
  - `pathband_post_xy`
  - `pathband_width_deg`
  - event frames (`occlusion_start_frame`, `occlusion_complete_frame`, `occlusion_end_frame`, `occlusion_end_complete_frame`)

This means the expensive geometric *planning* is not done during presentation.

### Runtime (per trial)
From `MoveDot1_experiment_occlusion_v14_PathBandOccluder.m`:
- Converts stored polylines from deg to px once per trial (`local_polyline_to_px`).
- Resolves draw options once per trial (`local_make_pathband_draw_options`).

### Runtime (per frame)
In the frame loop:
- Draws dot, arena, fixation.
- Draws pre/post occluder band with `local_draw_polyline_segments(...)` when active.
- In straight style, this function still:
  - applies terminal shifts,
  - rebuilds polyline from paired coords,
  - estimates tangents,
  - builds a closed polygon,
  - issues `Screen('FillPoly',...)` (plus optional start-width boost polygon).

So yes: polygon end handling is computed at runtime, every active frame.

## Expected frame-rate load

### Current default configuration helps
- `pathBandEntranceStyle = 'straight'` (default): uses a continuous polygon path, which is much cheaper than many segment quads/caps in round style.
- Pathband polyline span is limited to the occlusion segment (around deviance-to-occlusion-end), not full trial path.
- `Priority(MaxPriority(window))` is enabled.

### Main load multipliers
1. `debugShowFrameInfoOverlay = true` (DrawFormattedText every frame).
2. `debugShowFullPathOverlay = true` (full trajectory redraw every frame).
3. Round terminal mode (`pathBandEntranceStyle='round'`) increases draw-call count.
4. Extra trigger activity with DataPixx calls in-frame.

### Practical interpretation
- With debug overlays OFF and straight mode, 120 Hz is likely sustainable on a normal MEG stimulus machine.
- With both debug overlays ON, dropped frames become more likely.

## Trigger timing precision assessment

## What is good
- Condition/occlusion/type1 catch frame triggers are emitted at explicit frame indices in the frame loop.
- Trial output stores `frame_flip_time` and per-frame durations (`frame_duration_ms`) for post-hoc checks.

## Precision caveats
1. Trigger emission call (`local_emit_trigger`) happens before `Screen('Flip')` in each relevant frame.
   - This can introduce a lead/offset between TTL onset and actual photon onset.
2. No explicit logging of PTB `Missed`/`StimulusOnsetTime` return values from `Screen('Flip')`.
   - So dropped/deadline-missed frames are not explicitly flagged online.
3. Catch-question triggers are emitted with `frame=0` (time-based UI phase, not frame-locked trial motion).
4. `SkipSyncTests=1` is set unconditionally, reducing confidence in sync diagnostics.

## MEG interpretation
- If frame pacing is stable, trigger-to-event precision is typically within sub-frame regime, but still should be treated as requiring empirical validation (photodiode + trigger channel).
- For MEG analysis, this is usually acceptable only after confirming low dropped-frame rate and stable visual-trigger offset in the lab setup.

## High-priority TODO for next versions
- [x] 1. **Disable forced `SkipSyncTests=1` in MEG runs**.
   - Keep it only for explicit debug/development mode.
- [x] 2. **Add explicit flip diagnostics**.
   - Capture full `Screen('Flip')` returns (`VBLTimestamp`, `StimulusOnsetTime`, `FlipTimestamp`, `Missed`) and save per frame.
- [x] 3. **Add online dropped-frame guardrails**.
   - Warn/abort if missed flips exceed threshold. (warn is default, abort optional)
- [x] 4. **Precompute styled occluder geometry per trial**.
   - Move start/end terminal transforms and polygon assembly out of per-frame draw path when possible. Move it to before trial presentation starts.
5. **Harden trigger-to-visual alignment**.
   - Evaluate DataPixx video-synced write strategy (or trigger immediately after flip with measured timestamp and documented constant offset).
6. **Enforce MEG-safe runtime flags**.
   - Require explicit override with warning to enable `debugShowFrameInfoOverlay` and `debugShowFullPathOverlay` when `Conf.MEG=1`.
7. **Separate frame-locked vs non-frame-locked trigger classes in docs/output**.
   - Keep catch-question triggers explicitly marked as non-frame trial events.
8. **Add a short dedicated timing validation mode**.
   - N trials stress test, produce timing report (refresh stats, missed frames, trigger counts, worst-case frame time).

## Medium-priority TODO
1. Cache per-trial pre/post occluder draw artifacts (e.g., polygon vertices in px) to reduce repeated memory allocation.
2. Reduce trigger call overhead in high-event sections (minimize repeated DataPixx schedule reconfiguration where possible).
3. Add an automatic summary at end of run: median/95th/max frame duration and count of frames above 1.0*ifi and 1.5*ifi.

## Suggested MEG acceptance criteria (practical)
1. Dropped/missed flips: near-zero (ideally 0) in production run length.
2. Visual-trigger offset (photodiode vs trigger channel): stable and documented; jitter much smaller than one frame.
3. Occlusion event frame integrity: `occlusion_complete_frame` and `occlusion_end_frame` visually confirmed on representative trials.
4. Debug overlays disabled in all recorded sessions.

