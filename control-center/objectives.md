# 0. Define experiment
- [ ] Hyphothesis:
  - [ ] Step 1. Formulate questions & first experiment ideas
  - [ ] Step 2. Literature review on the specific question's subject
  - [ ] Step 3. Refine questions in light of literature
- [ ] Design & methods
  - [ ] Step 1. Decide a design and analysis to investigate questions
  - [ ] Step 2. Define specific hyphothesis and operationalize them into statistically testable outcomes
  - [ ] Step 3. Assess feasibility of specific design choices: i.e. can stimuli be generated in such a way that dRSA can work and allow testing the hyphotheis? -> model/models dRSA autocorrelation and cross-correlation have a structure that won't introduce biases or render hyphothesis untestable.  
  - [ ] Step 4. Literature review on the specific design & methods (paradigm, stimuli, analysis type, imaging technique or behavioural task, sample size, ...) aimed at finding tweaks for the design that will likely reduce β (Type II error) by increasing the SNR.
  - [ ] Step 5. Refine design in light of literature review

# 1. Prepare a presentation for MEG meeting 04/02/2026 - target audience: strong knowledge in cognitive neuroscience and research (PhD level at least up to full professor), good knowledge of EEG and/or MEG, no knowledge of dRSA but some knowledge of RSA expected.
- [ ] Prepare brief intro slides to frame in the literature and ends naturally flowing into research questions
- [ ] - [ ] prepare script for producing stimuli example for condition and figures of expected results for every hyphothesis and save as images
- [ ] a slide on the design with image of stimuli
- [ ] a slide showcasing the timing and structure of paradigm with trial with various presented screens
- [ ] Slides on: hyphothesis and expected result
- [ ] Slide showing some of the remaining open questions on the design (choose the ones that the specific audience could be able answer: again, audience knowladgeble in neuroscience and EEG/MEG but not in dRSA)

# 2. Implement final version of experiment
- [ ] define final stimuli generation code and parameters - e.g. of checks: dRSA autocorr and cross-corr
- [ ] define conditions and trial structure (runs, blocks, etc) - e.g. of checks: total exp time
- [ ] define experiment code - e.g. of checks: triggers sent correctly and carry info necessary for analysis

# 3. Run experiment pilot
- [ ] Run a single participant
- [ ] Prepare and test minimal analysis pipeline
- [ ] Run pilot of few subjects
- [ ] Prepare and test full analsys pipeline
- [ ] Implement adjustments if needed
  
# 4. Run full experiment
- [ ] Run rest of participants




## Deep Research

```text
You are conducting a theory-validity and hypothesis-refinement audit for a near-locked MEG experiment.
Your job is to determine whether the core scientific question and hypotheses are conceptually valid, malformed, or naive given current literature, then rewrite them into sharper, testable versions that fit the existing design and analysis.

Non-negotiable scope:
- Do NOT redesign the paradigm.
- Do NOT propose hypotheses that require new tasks, new modalities, or major protocol changes.
- Treat this as a single-moving-dot core question for inference (analyze one dot stream as primary; treat two-dot attention effects as optional secondary hypotheses only).

Project-specific design to anchor to:
- Modality: MEG.
- Stimulus: smooth moving dot trajectories in a bounded 2D rectangle (screen).
- Trial duration: 2.67 s.
- Sampling/display: 120 Hz (~320 frames/trial).
- Core conditions:
  - Non-deviant: smooth trajectory with constant curvature, no abrupt deviant turn.
  - Deviant: abrupt heading change at fixed midpoint onset (~1.335 s) and deviant-only post-onset curvature modulation.
- Predicted baseline for deviant trials is available as a no-deviant counterpart path (same pre-deviant-onset trajectory), enabling predicted-vs-observed contrasts after deviant onset.
- Analysis framework is dRSA-based (see "de Vries & Wurm, 2023.pdf" in the project for the dRSA framework) and mostly fixed.

Terminology convention (must use unless revised by evidence):
- Use "anticipatory representation" (not "predictive lag").
- Use "sensory-following representation" (not "lagged representation").
- Use "temporal offset" as the metric name, with sign convention:
  `Delta_t = t_brain - t_stim` (negative = anticipatory, positive = sensory-following).
- If literature indicates these labels are misleading or suboptimal, propose corrected terminology and apply it consistently.

Available analysis objects/readouts (stay within these):
- Non-deviant trials:
  - Primary model: observed-path position.
  - Secondary model (only if literature + separability support it): observed-path direction.
- Deviant trials:
  - Primary models: predicted-path position (no-deviant baseline), observed deviant-path position, PE-like position model (predicted minus observed).
  - Secondary models (only if separable): analogous direction models.
- Main outcomes: dRSA temporal-offset peak (`Delta_t_peak`) and dRSA peak amplitude over time.
- Existing constraints:
  - Position-vs-direction collinearity is currently problematic, so primary inference should rely on position models.
  - Autocorrelation limits temporal precision.
  - PE-model autocorrelation is currently high, so sharp timing claims for PE require caution.

Current hypothesis set to audit and rewrite:

Important framing constraints for all hypotheses:
- All hypotheses must be testable with current dRSA outputs (`Delta_t_peak`, peak amplitude, and their evolution across time windows).
- Separate claims about temporal-offset shifts from claims about representation-strength changes.
- When possible, use matched pseudo-onset in non-deviant trials for post-onset deviant contrasts.

H1a (non-deviant): within-trial emergence of anticipatory representation
- In predictable non-deviant trials, representation evolves from clearly sensory-following early in trial toward less sensory-following and potentially anticipatory later in the same non-deviant trials.
- Illustrative pattern to audit: early best match around clearly positive offsets (e.g., around +100 ms), later shift toward near-zero or negative offsets.
- Core interpretation under test: as trajectory parameters become more certain within trial, internal prediction strengthens.

Competing H1b (non-deviant): no anticipatory emergence
- Representation remains predominantly sensory-following and approximately stable across time.
- Any apparent offset drift is small/non-systematic or better explained by temporal-smearing/autocorrelation limits.

H2a (non-deviant): dampening of sensory-following representation strength
- As anticipatory representation strengthens, sensory-following peak amplitude decreases over time (explained-away account).

Competing H2b (non-deviant): sharpening of sensory-following representation strength
- Sensory-following peak amplitude increases over time (sensory-gain/tuning account).
- This is intentionally orthogonal to H1: there might be multiple peaks (i.g. an anticipatory and sensory-following one).

H3a (deviant, post-onset): PE-to-update sequence
- Time-lock to deviant onset.
- Early post-onset window: PE-like model transiently dominates over both predicted-path and observed-deviant-path models.
- Later post-onset window: representation shifts toward observed deviant-path model (update to new trajectory).

Competing H3b (deviant, post-onset): persistence-then-switch without strong PE dominance
- Early post-onset period is dominated by persistence of pre-deviant predicted-path representation.
- Transition to deviant-path representation occurs later, with weak or absent transient PE dominance.

Optional H4 (only if stable-vs-volatile context is actually implemented): context modulation of H1
- In non-deviant trials, stable context should show stronger/faster H1 emergence than volatile mixed context.
- Competing view: no context modulation or opposite modulation.

Optional H5 (secondary only): attention modulation
- If a second dot is presented (design choice), attended stream should show earlier and/or stronger representation than unattended stream.
- Treat as secondary due current single-dot core inference priority and experiment design; adding dot 2 is just an option in case there is not enough meat with other hyphothesis or as following experiment.

Primary goals:
1. Validate whether each hypothesis is fundamentally supported by literature assumptions.
2. Identify malformed or weakly grounded hypotheses (non-falsifiable, circular, over-claimed, or mismatched to measurements).
3. Rewrite hypotheses into precise, falsifiable, analysis-compatible statements.
4. Define boundary conditions and interpretation limits so claims stay conservative and defensible.
5. Recommend only minimal wording/contrast refinements compatible with current implementation.
6. Identify a minimal core hypothesis set that could support a strong first publication on its own.
7. Audit and, if needed, correct terminology so constructs are theory-valid, non-circular, and consistently mapped to measurable quantities.

Required output format:
A) Verdict summary (max 250 words): Is the core question/hypothesis logic sound overall?
B) Assumption audit table:
   Columns = Hypothesis | Core assumptions | Supporting evidence strength | Contradictory evidence | Verdict (Sound / Needs revision / Malformed).
C) Malformed-or-naive diagnosis:
   Explicit list of problems in current wording and logic.
D) Keep/Rewrite/Drop table:
   One row per current hypothesis item (H1a/H1b/H2a/H2b/H3a/H3b/H4/H5) with a short reason.
E) Rewritten prereg-ready hypotheses:
   For each rewritten H, provide:
   - exact prediction,
   - exact competing alternative,
   - exact falsification criterion,
   - minimal claim allowed if supported.
F) Analysis-compatibility check:
   For each rewritten H, state whether it is testable with current dRSA outputs and what exact contrast/readout is needed.
G) Evidence-grounded timing guidance:
   Literature-supported expected windows/signatures (if uncertain, say underdetermined).
H) Terminology audit and harmonization:
   - Proposed final term set
   - Definitions
   - Mapping from old -> new labels
   - Rationale from literature
   - Any residual ambiguity
I) Annotated references (20-40 key sources):
   Primary papers + strong reviews/meta-analyses; include DOI/link and 1-2 lines on relevance.

Quality rules:
- Include both supporting and conflicting/null results.
- Distinguish evidence-backed conclusions from inference/speculation.
- If evidence is insufficient, explicitly state "currently underdetermined".
- Avoid generic predictive-coding rhetoric; tie every claim to the literature and this exact design and measurable contrasts.
```
