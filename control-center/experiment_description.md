# NOTE (2026-02-24)
This file contains early brainstorming notes and open questions.

For the integrated, work-oriented version of this document (hypotheses,
design constraints, implementation map, and pre-data quality gates), start from:
- `control-center/experiment_blueprint.md`

# Theoric aim

Most of the neuroscientific investigation has been performed on responses to static stimuli. Not much is know about how the brain represents dynamic stimuli. Sometimes, when dynamic stimuli has been used, MVPA has been used to assess the decodability of different properties of the stimuli. Here we are revolutionasing the field by directly investigating the way dynamic stimuli are rapresented and how these representations evolve over time. Moreover, we can address questions linked to PC, in contexts that are more akin to the real computations our brain needs to do. If we walk down a street or we drive we need to continuosly represent and likely predict the movement of the objects that surround us; consequently the representation in our brain will continuosly evolve with new incoming information. This work will use Magnetoencephalography to record brain activity during the presenation of dots moving on a screen. Modelling togehter the stimulus and the brain activity we will look at when and which representations of the stimuli can be traced in the brain activity. Does the brain predict, or observe moving objects? Furthermore, manipulating the predictability of the stimuli with a sudden change in curvature (deviance), we will be able to directly test core hyphothesis of the Predictive Coding (PC) Theory  under circumstances that reflects the constantly evolving world around us, where unpredictable events determines changes in the ongoing dynamic world around us and force us to reshape our world model in real time. Using dRSA we can match models of the stimulus and brain data and the evolution of the peak of this match (peak strenght of dRSA), so we can look and what happens to the representation of the stimulus in the brain: can we see any predictive representation (brain represents the stimulus before it happens)? or we see lagged representations (brain representation follows the stimulus)? How does these representations evolve in time? from the appearance of the dots on, does the lagged representation get stronger or weaker? PC predicts that with time passing (during the trial) and the processing of path curvature and speed, the prior should strenghten, and consequently lagged representation should become weaker (dampening), but the sharpening hyphothesis would predict the lagged representation becomes stronger becasue it is sharpened. Another hyphothesis can be tested when the deviance happens (abrupt instantaneous change in direction and curvature): according to PC the brain should represent the prediction error (PE) until that causes the brain to re-adjust to a new prediction (the new direction/curvature); or does the brain continue to represent the predicted path and then switches to the new one (deviant), once it can realibly represent the new velocity and curvature? With dRSA we can model the two paths (predicted, deviant) and also model the PE, and test the timing of this dynamic evolution of brain representations. 


# How to build the actual stimuli to achieve the Theoric aim
 Stimuli generation is in place following some principles. All principles can be challenged and changed if there is a theoretical or practical reason to do it.

**LEGEND: Design choice --> reason**

- every trial will show a path (dot moving on the screen)
- there will be 20 unique paths per condition (deviant/non-deviant) repeated many times (20) randomly over the course of the experiment (400 trials per condition) --> to balance the need of variety in stimuli for generalization and avoid learning with the need of repetitions to reduce the noise in representations.
- the path will be of constant curvature and speed --> this allows for the formation of a strong prior/prediction/representation on how the path will evolve within the first few hundreds of millisenconds
- the 20 unique paths have all the same curvature? --> YES: makes formation of strong prior/inference of trajectory parameters quicker, stronger dRSA peaks sooner; NO: (1) more variability in the trials for diminishing learning, (2) makes slower the inference of trajectory parameters, so representation evolution might be seen better, but if too slow, it might impeed reaching of a strong enough peak? The aim is to reach a strong peak before the deviance occurs (right now at deviant time 1.33 middle of the trial), so it might be enough time even with curvature and speed change at every trial.
- At deviant point change curvature or just the instantaneous direction? We can maintain the same free parameter at the start of the trial and at the start of the deviant --> for a general sense of consistency no real spelled out reason. 
- speed? I would not mess up with it --> becasue I feel it changes also the representation correlation structure a lot? **Specific reason to be better defined**
- Length (in seconds): 2.66 s, speed: , --> it feels long enough and not too slow or quick
- 

 - BONUS: two dots might be presented in every trial instead of one, mainly to add a potential question on attention and representation.
 - the two dots must not be too close to each other --> **Specific reason to be better defined**
 - 


# what do we need to do to get at the questions we want to test?
-First, we try to optimize with one dot, then we add the second one if feasible.
DOT1:
- Define the models: position, direction, PE models? others?
- Position model is easy and defined as x,y of the dot at every timepoint
- Position model autocorrelation (correlate with itself) well? Yes and no, yes because only a central diagonal of autocorrelation, but a lot spreaded, so high correlations for a too big chunk of off-diagonal values more than 0.5s in each direction to drop the autocorrelation below 0.2. Ideally we would like to reduce it as it is used to define the autocorrelation boarder and then regress the autocorrelation out and if the boarders are too far, then the eventual dRSA correlation between brain and position model will be more spreaded and timing precision lower!
- can we use direction with the current stimuli generation? The issue is the cross-correlation between position and direction. This prevents running effective PCR to regreess out the direction model when analysis dRSA between brain and postion, and to regreess out the position model when analysis dRSA between brain and direction. This PCR step is essential to claim indipendence of results, namely to claim the position model and brain dRSA correlation means a brain representation only of the position and not of the direction, and viceversa. Question: Can we optimize the stimuli generation so that the cross-correlation fades or we would need to change things that than breaks other aspects of the design?
- 
