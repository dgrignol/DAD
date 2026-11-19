# MEG protocol CIMeC Ingmar de Vries

-> a double ^^ preceeding text indicate a part in which the text needs updating.
-> The current v25 runtime uses scheduled trial onset, an all-trial onset guard, an immediate dot-offset baseline, and the photodiode. Confirm the physical sensor
   placement against the configured patch before every acquisition.

## Before subject

- Turn on Polhemus in prep room (minimum 15 minutes before using).
  - Metal off table.
  - Twist dial on speakers to turn on.
  - double check Polhemus is ON.
- Acquisition PC:
  - Pword: m3g.
  - Click `acquisition` icon on desktop if software not already open.
  - Under `Settings` and next to `Project` click `Change` and select `ActionPrediction`. For `Subject` select your ethical code e.g. `2019-022`. Then:
    - File > Load settings, select ^^`Name of saved settings to insert` and click `Apply`.
    - 1000Hz, ^^N channels (check photodiode in misc is present, check eyetracker is present as ^^tyepToInsert), `upright` position.
    - IMPORTANT: Check online high pass filter to only DC offset or off or 0.03 Hz.
  - Press GO! and stop again.
  - Tune (once a day):
    - Tools > Tuner.
    - Load previous parameters: most recent or similar helium level. If bad pick other previous file.
    - `Measure noise`. Should be 0-4.
    - If necessary click `Tune` and wait till noise levels look better (2-3 min). If a channel is below 0 press ctrl then click on it and heat it with Commands > Heat channel. Then stop the tuning. `Stop measurement` > `Stop collector`.
    - Save tuning with convention: YYYYMMDDup53.
    - Then File > Exit.
  - Look at channels with empty room, if noisy:
    - Change the scale - typical is 440fT.
    - Select it. Commands > heat sensor (twice) > reset channels.
  - check camera in mangetic shielded room and turn it on if off.
  - If scanner is acting weird (e.g. dialogue boxes with error messages popping up):
    - Restart system: `RestartAcqPrograms` icon on desktop. When prompted on the terminal window type confirm with `y`; when asked if reboot or just restart servers, choose servers with `s`. If this fails you have to manually reset the switch for the real-time computers (SCC) in the Electronics cabinet (talk to technician before touching the cabinets).

- Prep room:
  - Start preparation software (follow instructions on desktop. Pword: m3g).
  - Load Project and settings as on acq PC.
  - Check 6 clean electrodes and stickers, and headcoils.
  - Check consent forms (MEG, privacy, sensitive data), payment and MRI release forms, and the task instructions.

- mangetic shielded room:
  - Projector on.
  - Lights off.
  - Reverse screen: VPutil icon on desktop. Type `4` and hit enter, then `pm r` and hit enter.
  - Check if button press works.
  - Check eyetracker in correct place.
  - Battery for photodiode and turn box on.
  - Attach photodiode to screen.

- Stimulation PC:
  - Open Matlab.
  - Confirm the launcher uses
    `MoveDot1_experiment_occlusion_v25_dotOffsetBlank.m` with
    `lib/Config_runtime_v25_dotOffsetBlank.m`.
  - Before launching the experiment, run `display_refresh_rate_utility.m`
    separately and confirm that the projector/DataPixx path is configured for
    `120 Hz`. The experiment only measures and validates refresh timing; it
    does not configure the display.
  - At startup, confirm that the reported Psychtoolbox and DataPixx refresh
    estimates are both within `0.5 Hz` of the configured `120 Hz`. The v25
    runtime aborts outside this tolerance by default.
  - Confirm the startup log reports `SkipSyncTests=0` and the applied lab
    `SyncTestSettings`: maximum standard deviation `0.0005 s`, 50 samples,
    maximum deviation `0.1`, and maximum duration `5 s`. These settings are
    already applied by the lab runtime preset; do not paste separate PTB
    startup commands or add `Screen('CloseAll')` to the launcher.
  - Check the configured photodiode patch before attaching the sensor:
    `photodiodePatchSizePx` defaults to the physically validated `[28 25]`
    and `photodiodePatchPositionPx` to `[193 0]` from the PTB window top-left.
    Startup must report `rect=[193 0 221 25]`; place the sensor over that patch.
  - Start a test and verify both the digital triggers and the photodiode
    pulses on the acquisition PC (photodiode channel: MISC 008).
  - Verify that the patch returns to black after each pulse and at experiment
    cleanup. Defaults are four frames for trial-condition markers and two
    frames for marked visual events.
  - In a test recording, verify the one-time trigger-`153` warm-up at the
    beginning of the MATLAB experiment launch: 10 four-frame high-patch
    cycles, each followed by an explicit low frame (about `0.42 s` at
    `120 Hz`). No moving dot or task event should appear during this warm-up.
    Trigger `153` is a warm-up marker, not a trial-condition marker.
  - The completed v24 physical audit measured 433 condition onsets at
    `18.76 +/- 0.47 ms`, with no guard exhaustion, and passed stop/resume.
    V25 retains that onset path. Before participant acquisition, run a short
    visual check with the real EyeLink: the dot must disappear immediately
    after its final motion frame, and each completed run must print a
    `Dot-offset summary` near one frame with `positiveMissed=0`.
  - Inspect `trial_onset_guard_*` for every retained trial. A retry produces
    only an extra low-patch fixation frame and no condition trigger. Up to six
    retries are allowed. If exhaustion occurs, v25 continues and sets
    `trial_onset_guard_exhausted` plus
    `trial_onset_guard_continued_after_exhaustion`; inspect that event's
    photodiode latency individually. Zero exhausted guards is expected and is
    not a test failure.

- Eyetracker:
  - Turn on power supply (power plug next to mangetic shielded room).
  - Turn on PC and screen.
  - Check if works.

## Prepare subject
- Explain preparation, scanning and experiment to participant.
- They sign the consent forms and MRI release form.
- Double check:
  - Normal vision (or have contact lenses / know their prescription for the MEG glasses).
  - No metal in body e.g. dental wires, surgical pins/plates.
  - Did not go in MRI scanner today.
- Shoes and pyjamas and plastic socks in shoes and gloves and mask without metal.
- Remove all piercings, jewellery, bra, hairbands and clips etc.
- Remove makeup.
- Sit in chair and read task instructions.
- Participant number: DOB in YYYYMMDD, first and third letter of mother's first name, and first and third letter of mother's maiden name. e.g. 19900905LNDR.
- Head coils with the white tape in positions and colour coding indicated on the polystyrene head.
  - Flat side on the skin.
- EOG + ECG (if you use them).
  - EOG: above and below the left eye, centre of electrode roughly in line with their pupil when looking forward. One electrode each side outside eyes near the temples.
- Make sure project ^^`projectName`, protocol ^^(yourProtocol) and settings (yourSetting) loaded on Subject preparation software.
- Put the glasses on the subject and take care not to move them while you digitize points.
- Click `Change` next to HPI, then `Coordinate frame alignment`.
- Digitize anatomical landmarks: the nasion and the two pre-auricular points.
- Check landmarks: hold pen away from head and press button, then on PC click `Coordinate frame alignment` again and select `Check`, take the points again - they should be almost the same as your first go. If big difference `Coordinate frame alignment` > `Redo`.
- Digitize head coils in order indicated on polystyrene head their left to right.
- Digitize around 250-300 points to give the shape of the head.
- Hold pen far from participant and click button till you hear a double beep.
- Click `ok`. Then File > Save preparation and make a note of the preparation file name.
- On acquisition PC, File > `Load normal preparation` and select the current file name.
- Gently remove the preparation glasses and gather up all the wires, being careful not to pull on them. Give them to the participant to hold.

## Recording
- Turn off fan (console is in side room on the right).
- Plug in EOG, ECG and head coils.
- Check photodidoe position and patch correspondance (util script present in MEG computer in Public folder), adjust patch positon if needed.
- Last explain to participant:
  - Relax in chair, not sit up too straight, find comfy position.
  - Show them button boxes: red = trajectory did NOT change, blue = trajectory changed.
  - Move chair up so they are touching the top of the helmet.
  - Ask them to take a big slow deep breath in and then out, they will slouch a little and then you put them up a tiny bit more in the chair so they are really up in the helmet.
  - Eye-tracker calibration.
  - Intercom.
  - Door opening in case of emergency.
  - Last short explanation important aspects task, including: relax, no eye movements, blink during blink breaks, fixate center.
- Check if eye-tracker sees eyes properly (start experiment, in eyetracker screen press ENTER and left or right arrow to see both eyes).
- On the acquisition PC, press Go! and check how the signal is. Type in 11 on intercom and ask them to blink a few times - check EOG and ECG channels if you have them, otherwise, look at eyetracker channel in MEG and frontal sensors (e.g. 923-922-921) in MEG. Ask them to take a big breath - should barely show up in the signal. Press Stop to end the recording.
- On stim PC, press C to calibrate eyetracker, then V to validate. Then `Accept` on eyetracker PC. Tell participant we are nearly ready to begin. Put intercom in only listening (^^how is it done in the new intercom?)
- On acquisition PC, click Go!
- Click yes to HPI measurement, if suggested to accept click ok.
  - Third number of head origin should be 40-50mm. Greater than 60mm is bad. **Note down the number**, the aim is for it to remain constant accross runs. At every run, that window will pop-up again and you can check the number did not change much.
  - If having problems with coils go in the mangetic shielded room and make sure they are attached properly, and that they are in the helmet, that you loaded the correct preparation.
- Tick `Record raw` before starting calibration.
- ToBeConsidered: Ask the participant to blink a few times and make horizontal and vertical eye movements before starting the experiment (it might become useful in the preprocessing for the ICA to better pick up components.)
- Keep acquisition and raw recording active before confirming the experiment
  start on the stimulation PC. This ensures that the one-time trigger-`153`
  warm-up and its photodiode pulses are recorded before the first participant
  trial.
- During recording, keep an eye on:
  - MEG signal, especially in your sensors of interest.
  - Note down noisy sensors (can heat them if necessary before next block).
  - Triggers and photodiode channel (MISC 008), EOC and ECG signals. Expect
    the one-time `153` warm-up sequence, followed by trial-condition triggers
    `71`, `81`, or `91`; do not score `153` as a trial.
  - Eyetracker monitor to check fixation.
- At the end of the block click Stop on acq PC and save the file using the proper naming system:
  - Date of birth, 1st and 3rd letter mum's first name and last name, current date plus time rounded to nearest 15 minutes, EC protocol (for us 2025059), run. For example:
    - 19890116CRRS_202003281125_2019022_run1
  - Blue changes; orange stays same.
  - Copy = CTRL + INS.
  - Paste = up arrow + INS.
- At the end of the block when you save, there is a pop-up window about the **head position** (to see if the head position changed), good values for the difference is anything <40/60mm (^^should be the same as before and difference smaller than mabe 20?)

## After recording
- Lower chair, unplug coils/electrodes. Take them to prep room and take off coils/electrodes. Let them change. Get them to fill out payment forms.
- Clean head coils and electrodes with alcohol and paper towels.
- Clean anything the subject came in contact with.

## End of day
- If misnamed files, rename like so:
  - Open a shell on the MEG console (right click in an empty space on desktop and select Konsole) and navigate to the local dataset repository for the day and your ethical code / project.
  - Linux command `ls` to show what is in the directory.
  - Now change each "fif" file name using the Linux command `mv oldFileName newFileName`.
- Turn off:
  - stim PC.
  - screen of acq PC.
  - camera monitor.
  - eyetracker PC.
  - photodiode box, remove battery and put in charger in prep room.
  - projector (first on remote then switch).
- Close mangetic shielded room door.
- Turn off Polhemus.
- Leave lab tidy.




## Troubleshooting
- If triggers or responses do not work in MEG lab:
  - Turn off everything and then turn on in correct order: pc, blue box, projector, flip projector with desktop icon and type `pm r`, then turn on eye tracker.
  - Note that the second screen will only turn on if the main blue box is on.
