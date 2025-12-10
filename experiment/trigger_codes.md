# Trigger Codes (MoveDot1_experiment_vX)

| Code / Range | Short label | Extended description |
| --- | --- | --- |
| 1–N (block 1, no-catch) | Trial onset block 1 | On first frame of a non-catch trial when attending `Conf.DotColor1`. Value is `condMatrixShuffled(2, iTrial)` (the `sequence` field from the stimulus file); if `condMatrixShuffled(1, iTrial) == 40` an extra +20 is added (legacy branch currently unused with the generated conditions). |
| 41–(N+40) (block 2, no-catch) | Trial onset block 2 | Same mapping as above but with +40 when `blockCondition == 2` (attend `Conf.DotColor2`). |
| 102 | Trial onset with catch | Overrides the onset code whenever the trial contains at least one catch; independent of block/condition. |
| 81 | Trial end (no-catch) | Sent on the last frame of trials without catches. |
| 103 | Trial end (with catch) | Sent on the last frame of trials that contain a catch. |
| 100 | Catch start (fix/occl) | First frame of a catch trial, used for both fixation and occlusion catches. |
| 101 | Catch end (fix/occl) | First frame after a catch finishes, used for both fixation and occlusion catches. |
| 113 | Missed fixation catch | First “Missed” frame when no response is given during a fixation catch; no equivalent code is emitted for occlusion misses (that branch is commented). |
| 150 | Gaze break | Emitted when the eye leaves the fixation window long enough to abort the trial; trial is queued for replay. |
| 151 | Replay start | Sent at the start of a replayed trial (after a gaze break) to mark the replay instance. |
| 201 | Response (left/“1”) | Button-response pulse `response+200` when the response value is 1; with current MEG button mapping all detected presses are recoded to 1, so 201 is the only response TTL actually sent. |
| 202 | Response (right/“2”, unreachable) | Would be emitted if `response == 2` ever occurred; current MEG code maps both button codes (1 and 8) to `response = 1`, so 202 is not produced in practice. |

**Notes**
- `N` depends on the `sequence` values in `condMatrixShuffled` (from `xySeqs(...).sequence`); with the current stimuli generator this is the per-condition trial index and typically ranges 1..trialsPerCondition (before the optional offsets above).
- All onset codes switch to 102 for trials containing catches, so the block/condition mapping only applies to pure non-catch trials.
