# Trigger Codes (MoveDot1_experiment_vX)

| Code / Range | Short label | Extended description |
| --- | --- | --- |
| 1–20 or 21–40 (blockCondition 1, no-catch) | Trial onset (attend DotColor1) | First frame of a non-catch trial when `blockCondition == 1`. Start from `sequence` (1–20). If `condition == 0` (directionVariance = 0; no deviation / more predictable), keep 1–20. If `condition == 45` (directionVariance = 45; deviation / less predictable), add +20 → 21–40. |
| 41–60 or 61–80 (blockCondition 2, no-catch) | Trial onset (attend DotColor2) | Same mapping as above, plus +40 for `blockCondition == 2`. So 41–60 for condition 0, or 61–80 for condition 45. |
| 102 | Trial onset with catch | Overrides the onset code whenever the trial has at least one catch; independent of block/condition. |
| 81 | Trial end (no-catch) | Sent on the last frame of trials without catches. |
| 103 | Trial end (with catch) | Sent on the last frame of trials that contain a catch. |
| 100 | Catch start (fix/occl) | First frame of glitching. |
| 101 | Catch end (fix/occl) | First frame after glitch finishes. |
| 113 | Missed fixation catch | First “Missed” frame when no response is given during a fixation catch; no equivalent code is emitted for occlusion misses (that branch is commented). |
| 150 | Gaze break | Emitted when the eye leaves the fixation window long enough to abort the trial; trial is queued for replay. |
| 151 | Replay start | Sent once, immediately before the first replayed trial begins (after a gaze break). A 100 ms pause follows the pulse before the replay trial proceeds. |
| 201 | Response (value 1) | Button-response pulse `response+200` when the response value is 1; the current MEG mapping recodes button codes 1 and 8 to `response = 1`, so 201 is the typical response TTL. |
| 202 | Response (value 2, NOT expected) | Would be emitted if `response == 2`; note that 202 is not in `TriggerValues`, so it will be flagged if it occurs (e.g., from a nonstandard button code). |

**Notes**
- `N` depends on the `sequence` values in `condMatrixShuffled` (from `xySeqs(...).sequence`); with the current stimuli generator `sequence` is 1–20, yielding the concrete ranges above.
- In the current stimuli generator, `condition` corresponds to `directionVariance` (predictability); Config.likelihood uses `[0, 45]`, where 0 means no deviant angle change at the deviant onset and 45 introduces a deviation.
- `blockCondition` is read from `BlockOrder` in the trial struct, so the attend-dot mapping can swap between blocks.
- All onset codes switch to 102 for trials containing catches, so the block/condition mapping only applies to pure non-catch trials.
