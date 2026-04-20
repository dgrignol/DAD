TODO
Current set of experiment script and config:
experiment/lib/Config_stimuli_generation_V28_rescueTraject.m
experiment/stimuli_generation_V28_rescueTraject.m
experiment/lib/Config_schedule_CreateInputV21_MoveDotV18_rescueTraject.m
experiment/CreateInputFiles_v21_rescueTraject.m
experiment/lib/Config_runtime_v18_rescueTraject.m
experiment/MoveDot1_experiment_occlusion_v18_rescueTraject.m




**Experiment Flow**
- [ ] accuracy printed at the end of every block
- [ ] Debug the current startup EyeLink calibration, which does not seem to work reliably.
  - [ ] at every re-calibration at the end of a block, it stops the experiment because it requires to press esc at the end of the calibration, this does not happen at the end of the very first calibration when you start the exp.

**From participant feedback**
- [ ] consider insertion of higher number of catch trials

**DONE**
**Procedure**
- [x] make PRACTICE! a non meg script that launches equal for everyone, only catch trials, tunable number of trials per run1 and 2/3, clear separation between run 1 and 2/3 with instructions of one and the other; Repeatable with a button press; it just ask the particpant number at the start to save a single file with the accuracy and RTs; No MEG part, but there is the eyelink part with calibration and replays but can be switched off with parameter at the start.

**Experiment Flow**
- [x] Add an option to call EyeLink calibration again between blocks.
- [x] Confirm that replayed trials should stay within the same run. The current behavior is run-local: replays are appended to the end of the current run (`run 1 = always_visible`, `runs 2/3 = occluded`). -> They do replay per run. Replay scheduling is run-local: when a fixation break (trigger 150) occurs, that trial’s source index is appended to the end of the current run’s schedule, and replay state resets at each run (no cross-run/block carry-over). With v19-generated TrialOrder, run 1 sources are always_visible and runs 2/3 sources are occluded, so replayed trials inherit the current run’s condition family; this condition split depends on the input schedule.
- [x] Decide the EyeLink calibration policy. The current v16 runtime calls calibration once at experiment start, before the block loop.
- [x] Add a visible countdown to the run-1 to run-2 transition screen (`End of first task. You can rest your eyes. Second task will start soon.`).
- [x] Add a trigger for experiment termination via `Esc`.
- [x] Standardize the `Esc` abort messages so they include block number and context, at minimum:
      - `esc pressed during block [block_number]`
      - `esc pressed at the end of block [block_number]`
- [x] Keep the current overwrite refusal, but add the help message: `if you intended to run from a specific block, set iBlock variable`
- [x] Review and finalize post-deviance trajectory controls in `Config_stimuli_generation_V26_eyeTrackerReplay.m`: `initialCurvatureWindows`, `deviantCurvatureWindows`, `deviantSignedTurnWindows`, `likelihood.directionChange`, `flipCurvatureOnDeviant`, and `randomizeCurvatureOnDeviant`.

**Saving and Output**
- [x] Keep the current block-level `.mat` save, which already happens before the end-of-block screen.
- [x] Make the EyeLink output block-independent as well, so all files needed for analysis exist by the time the end-of-block message appears.
- [x] Decide whether EyeLink should use one EDF per block. The current code opens one EDF before the block loop and receives/moves it only once during final cleanup.

**Visual Display**
- [x] Set `background`, `rectColor`, and `rectBorderColor` to the same black color: `[0 0 0]`. The current mismatch is mainly the dark-gray arena fill against a black background.
- [x] Check the full-screen mask appearance after the color change. The post-trial grid mask still uses non-black colors. -> mask has been removed for now.
- [X] If black side bands are still visible, verify whether they come from display/fading rather than explicit band drawing in the v16 code. -> not needed if background is black [0 0 0]













OLD STUFF TO BE IGNORED FOR NOW:
- [x] check changes to constrains worked: new constrain that checks the non-deviant paths in deviant trials would not have touched the boudaries
- [x] scale the jitter with the distance from the center? (Bence suggestion)
- [ ] Catch trials behavioural testing
- [x] try to find if eye-tracer detects blinks automatically or not and adjust code so that it is less demanding on the participant
- [ ] non avere anticorrelation della posizione dei dots.

- [ ] SAVING ISSUE:
    - [ ] extend mid-run save "firstBlock" to the end too and call it "backupSave"
    - [ ] at the end see if you can save first and then send the message to finish the experiment
    - [ ] Check this issue carefully, especially the : CatchOutput construction can crash right before saving. CatchOutput is built before the final save(...) (experiment/MoveDot1_experiment_vX.m (lines 2413-2434)). It assumes that for each trial iCatch = 1:numcatch, arrays like correct_response and RT have at least numcatch elements (experiment/MoveDot1_experiment_vX.m (lines 2427-2428)). If a trial ends early (e.g., fixation break / replay logic) and you don’t end up logging responses for all planned catch trials, this can throw an indexing error at the very end → no full-run file. 
    - [ ] check also this potential issue: Teardown happens before saving. The full-run save is after EyeLink shutdown and DataPixx close (experiment/MoveDot1_experiment_vX.m (lines 2355-2383)). Those calls aren’t fully protected by try/catch (only ReceiveFile is), so any error there would skip the final save.
    - [ ] If esc pressed, should the experiment save a partialData file before exiting?
- [ ] update the scripts to the new folder structure with something like "I've changed the folder structure in data and I would like to update the scripts inspect_fif_report.py accordingly
If I ask for --subject 4 --run 1 the scripts should take
data/sub04/MEG/sub04_run01.fif"
- [ ] Create script for fixing the issue with triggers seen in the simulation with 99: a few expected triggers were sent just after an unexpected one that impeded the mne.find_events function to correctly parse the trigger.
- [ ] Check triggers of the real run. For sure there will be a problem with the currect script that does not take into account possibility of replays (so 150 instead of 81 to end a trial and then other unexpected trials at the end). Use the 150 and 151 for adjusting expectations and create a report section of replayed trials with some visualization and summary stats (number of replays, when in the experiment timeline, how many replays per gaze break angle setting)

DONE
FROM DAVIDE:

- [x] 1440 screen size (width). Metti bande nere per renderlo 1440 da 1960.
- [x] 1/conf.refrets

- [x] Try experiment practice. Some trials seems not to register the response. Issue of timing? If you respond too close to the end of the trial it does not collect the response? It feels to much now. Or is it designed on puprpouse like this to increase difficulty (i.e. to push RTs?) -> by design the end of the trial does not collect response (0.2 sec)
- [X] Implement new logic for detecting drifts or saccades away from the fixation cross, stop the current trial and repeat that or append it to a list of trials to re-run at the end of the block. (Second solution better but more difficult?). Implement a switch off for this.
- [X] feedback for false positive quando dice che  c'è ma non c'è
