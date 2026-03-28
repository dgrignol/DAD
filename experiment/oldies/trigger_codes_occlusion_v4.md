# Trigger Codes (MoveDot1_experiment_occlusion_v4_catchTrials)

| Code | Short label | Extended description |
| --- | --- | --- |
| 31 | Trial onset (always_visible) | First frame of an `always_visible` trial. |
| 41 | Trial onset (occluded_nondeviant) | First frame of an `occluded_nondeviant` trial. |
| 51 | Trial onset (occluded_deviant) | First frame of an `occluded_deviant` trial. |
| 111 | Occlusion start | First frame where the dot starts disappearing (`occlusion_start_frame`). |
| 112 | Occlusion complete | First frame where the dot is fully invisible (`occlusion_complete_frame`). |
| 114 | Occlusion end start | First frame where the dot starts reappearing (`occlusion_end_frame`). |
| 115 | Occlusion end complete | First frame where the dot is fully visible again (`occlusion_end_complete_frame`). |
| 116 | Catch type-1 disappear | Type-1 catch only: first frame where dot is hidden in run-1 disappear/reappear catch. |
| 117 | Catch type-1 reappear | Type-1 catch only: first frame where dot reappears after invisible interval. |
| 118 | Catch question start | First frame of the question prompt: `Has the dot changed its path?` |
| 101 | Catch question end | Emitted when catch question closes (response or timeout). |
| 201 | Catch response YES | Participant answered YES (`changed`). |
| 202 | Catch response NO | Participant answered NO (`not changed`). |
| 113 | Catch timeout | No response within configured catch response window. |

**Notes**
- This trigger map is specific to `experiment/MoveDot1_experiment_occlusion_v4_catchTrials.m`.
- Catch planning is generated upstream in `CreateInputFiles_v14_threeRunsPerBlock_catch.m` and consumed at runtime from `CatchPlan`.
- Type-1 catches are inserted in run 1 and suppress occluder rendering.
- Type-2 catches are inserted in runs 2 and 3 and keep normal occlusion rendering; only the end question is added.
- Response collection supports MEG button-box input and keyboard debug input.
