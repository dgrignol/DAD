Questions:
- [x] Try experiment practice. Some trials seems not to register the response. Issue of timing? If you respond too close to the end of the trial it does not collect the response? It feels to much now. Or is it designed on puprpouse like this to increase difficulty (i.e. to push RTs?) -> by design the end of the trial does not collect response (0.2 sec)
- [X] Implement new logic for detecting drifts or saccades away from the fixation cross, stop the current trial and repeat that or append it to a list of trials to re-run at the end of the block. (Second solution better but more difficult?). Implement a switch off for this.
- [ ] scale the jitter with the distance from the center? (Bence suggestion)
- [ ] Catch trials behavioural testing
- [X] feedback for false positive quando dice che  c'è ma non c'è
- [ ] try to find if eye-tracer detects blinks automatically or not and adjust code so that it is less demanding on the participant
- [ ] non avere anticorrelation della posizione dei dots.


FROM DAVIDE:

- [ ] 1440 screen size (width). Metti bande nere per renderlo 1440 da 1960.
- [ ] 1/conf.refrets