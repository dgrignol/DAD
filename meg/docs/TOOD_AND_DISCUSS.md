IMPORTANT point to add/think about

Analysis interpretability & documenting
->the first part of every html report should contain the info about the parameters and settings commented in a way that each one is understandable (how to do? is json already printed at the start? add comments and more info in json that is already printed at the start?)

PREPROCESSING:
->Notch filter? (present for Ingmar not in FLUX) -> solution cleanline style filter or lowpass at 40 - Simple: Sinusoidal regression / DFT filtering: Good when frequency and amplitude are reasonably stable. MNE offers method="spectrum_fit", set it for the line noise freq
    ->> Make both solutions as other options and save all the files (create tree folder structure to accomodate and then to use as input from dRSA and compare results)
->all data cleaning in separate step? (cleaner but feasible for us? it is just a matter of propagating the info that we should anyway doing, so maybe "06 cleaning and epoching"? I think we are already close to that, just the ICA components are removed before?)
->smoothing necessary? why? (the downsample is already "smoothing")
->NaN bad-chans convention better than removing them? (it is cleaner but it adds a step before inserting the data in the dRSA pipeline?)
-> Address deviation from flux "headshift", can we adjust the head over runs?



DECODING:
--> decoding of the occluded condition should see bump at reappearence
-> ERF of occluded condition
-> decoding of the position (x,y) itself not just left-right


DONE
->!!!!save eyetracking!!!!!!
-> ITI from [500,1000] to [1000,1500]? Ingmar suggestion. Plus, every 10 trials a mini break of 5 seconds to blink with a written text made for blinking?
