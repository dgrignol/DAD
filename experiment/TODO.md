Questions:
- [x] Try experiment practice. Some trials seems not to register the response. Issue of timing? If you respond too close to the end of the trial it does not collect the response? It feels to much now. Or is it designed on puprpouse like this to increase difficulty (i.e. to push RTs?) -> by design the end of the trial does not collect response (0.2 sec)
- [X] Implement new logic for detecting drifts or saccades away from the fixation cross, stop the current trial and repeat that or append it to a list of trials to re-run at the end of the block. (Second solution better but more difficult?). Implement a switch off for this.
- [ ] scale the jitter with the distance from the center? (Bence suggestion)
- [ ] Catch trials behavioural testing
- [X] feedback for false positive quando dice che  c'è ma non c'è
- [ ] try to find if eye-tracer detects blinks automatically or not and adjust code so that it is less demanding on the participant
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


FROM DAVIDE:

- [x] 1440 screen size (width). Metti bande nere per renderlo 1440 da 1960.
- [x] 1/conf.refrets