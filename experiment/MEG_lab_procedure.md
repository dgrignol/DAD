# MEG protocol CIMeC Ingmar de Vries

## Potential issues
- If triggers or responses do not work in MEG lab:
  - Turn off everything and then turn on in correct order: pc, blue box, projector, flip projector with desktop icon and type `pm r`, then turn on eye tracker.
  - Note that the second screen will only turn on if the main blue box is on.

## Before subject
- Turn on Polhemus in prep room (minimum 15 minutes before using).
  - Metal off table.
  - Twist dial on speakers to turn on.

- Acquisition PC:
  - Pword: m3g.
  - Click `acquisition` icon on desktop if software not already open.
  - Under `Settings` and next to `Project` click `Change` and select `ActionPrediction`. For `Subject` select our ethical code `2019-022`. Then:
    - File > Load settings, select `std` and click `Apply`.
    - 1000Hz, 336 channels (check photodiode in misc), `upright` position.
    - IMPORTANT: Change online high pass filter to off or 0.03 Hz or only DC offset.
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
  - Turn on camera in MSR.
  - If scanner is acting weird (e.g. dialogue boxes with error messages popping up):
    - Restart system: `RestartAcqPrograms` icon on desktop. When prompted on the terminal window type confirm with `y`; when asked if reboot or just restart servers, choose servers with `s`. If this fails you have to manually reset the switch for the real-time computers (SCC) in the Electronics cabinet (best to talk to Davide and Gianp before touching the cabinets).

- Prep room:
  - Start preparation software (follow instructions on desktop. Pword: m3g).
  - Load Project and settings as on acq PC.
  - Check 6 clean electrodes and stickers, and headcoils.
  - Check consent forms, payment and MRI release forms, and the task instructions.

- MSR:
  - Projector on.
  - Lights off.
  - Reverse screen: VPutil icon on desktop. Type `4` and hit enter, then `pm r` and hit enter.
  - Check if button press works.
  - Check eyetracker in correct place.
  - Battery for photodiode and turn box on.
  - Attach photodiode to screen.

- Stimulation PC:
  - Open Matlab 2012b.
  - Start experiment and check triggers and photodiode on acquisition PC.
  - Check and maybe change refresh rate of projector to 60 instead of 120.

- Eyetracker:
  - Turn on power supply (power plug next to MSR).
  - Turn on PC and screen.
  - Check if works.

## Prepare subject
- Put on white apron before subject arrives.
- Explain preparation, scanning and experiment to participant.
- They sign the consent form and MRI release form.
- Double check:
  - Normal vision (or have contacts / know their prescription for the MEG glasses).
  - No metal in body e.g. dental wires, surgical pins/plates.
  - Did not go in MRI scanner today.
- Shoes and pyjamas and plastic socks in shoes and gloves and mask without metal.
- Remove all piercings, jewellery, bra, hairbands and clips etc.
- Remove makeup.
- Sit in chair and read task instructions.
- Participant number: DOB in YYYYMMDD, first and third letter of mother's first name, and first and third letter of mother's maiden name. e.g. 19900905LNDR.
- Head coils with the white tape in positions and colour coding indicated on the polystyrene head.
  - Flat side on the skin.
- EOG + ECG.
  - EOG: above and below the left eye, centre of electrode roughly in line with their pupil when looking forward. One electrode each side outside eyes near the temples.
- Make sure project `OscillatingWM`, protocol (2009-036) and settings (std) loaded on Subject preparation software.
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
- Last explain to participant:
  - Relax in chair, not sit up too straight, find comfy position.
  - Show them button boxes: red = match, blue = no match.
  - Move chair up so they are touching the top of the helmet.
  - Ask them to take a big slow deep breath in and then out, they will slouch a little and then you put them up a tiny bit more in the chair so they are really up in the helmet.
  - Eye-tracker calibration.
  - Intercom.
  - Door opening in case of emergency.
  - Last short explanation important aspects task, including: relax, no eye movements during memory encoding, blink in between trials.
- Check if eye-tracker sees eyes properly (start experiment, in eyetracker screen press ENTER and left or right arrow to see both eyes).
- On the acquisition PC, press Go! and check how the signal is. Type in 11 on intercom and ask them to blink a few times - check EOG and ECG channels. Ask them to take a big breath - should barely show up in the signal. Press Stop to end the recording.
- On stim PC, press C to calibrate eyetracker, then V to validate. Then `Accept` on eyetracker PC. Tell participant we are nearly ready to begin then hold down `T` on intercom (now you hear them but they do not hear you. Press `T` twice for them to hear you again).
- On acquisition PC, click Go!
- Click yes to HPI measurement, if suggested to accept click ok.
  - Third number of head origin should be 40-50mm. Greater than 60mm is bad.
  - If having problems with coils go in the MSR and make sure they are attached properly, and that they are in the helmet, that you loaded the correct preparation.
- Tick `Record raw` before starting calibration.
- During recording, keep an eye on:
  - MEG signal, especially in your sensors of interest.
  - Note down noisy sensors (can heat them if necessary before next block).
  - Triggers and photodiode channel (MISC 008), EOC and ECG signals.
  - Eyetracker monitor to check fixation.
- At the end of the block click Stop on acq PC and save the file using the proper naming system:
  - Date of birth, 1st and 3rd letter mum's first name and last name, current date plus time rounded to nearest 15 minutes, EC protocol (for us 2009036), run. For example:
    - 19890116CRRS_202003281125_2019022_run1
  - Blue changes; orange stays same.
  - Copy = CTRL + INS.
  - Paste = up arrow + INS.

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
- Close MSR door.
- Turn off Polhemus.
- Leave lab tidy.
