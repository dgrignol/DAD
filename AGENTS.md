# AGENTS.md

## Cohort Execution Policy
- Run independent participant operations in parallel by default. Submit the
  whole eligible cohort together and use the available allocated CPU/memory
  resources; do not silently process a parallelizable cohort sequentially.
- Use separate participant processes or scheduler array tasks with disjoint
  outputs. Limit nested numerical-library threads to avoid oversubscription.
- Keep genuine dependencies and shared final publication ordered. In particular,
  complete whole-table choice validation before any bulk import writes, and join
  participant results before publishing cohort summaries or scientific reviews.
- Preserve manual work, source/model validation, review gates and immutable
  outputs. Scientific acceptance does not follow from parallel job completion.
- Expose worker/resource controls, document unavoidable serial steps, and print
  participant/run progress. This policy also applies to future cohort workflows.

## Documentation and Commenting Requirements
- Add extensive documentation at the start of every script with example of usage.
- When adding or modifying code, include section-level comments that explain the purpose of each major block and how data flows through it.
- Add concise inline comments for non-obvious logic, edge cases, or
  experiment-specific assumptions; avoid redundant "what" comments.
- Keep comments consistent with the existing file style and use ASCII unless
  the file already uses non-ASCII.
- Update or add docstrings for scripts/functions when behavior changes,
  covering inputs, outputs, and key assumptions.

## Local Tooling
- MATLAB is installed on this machine at:
  `/Applications/MATLAB_R2020a.app/bin/matlab`
- Always invoke MATLAB via the absolute binary path. Do not use `matlab`
  from `PATH` in commands.
- Example non-interactive test run (preferred):
  `/Applications/MATLAB_R2020a.app/bin/matlab -batch "run('experiment/stimuli_generation_v12.m')"`
- If host architecture detection fails on Apple Silicon, force Intel-arch
  startup for R2020a:
  `arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch "run('experiment/stimuli_generation_v12.m')"`

## Stage-07 Group Report Policy
- Group linked and embedded HTML reports must contain only group-average
  figures. Never embed participant galleries or participant Stage-08/09 plots
  in a group report. Keep individual figures in participant reports; JSON
  inventories may retain their provenance.
- Use portable relative image links and show a visible clickable figure path.
- For presentation-only maintenance, use Stage 07's
  `--group-html-only` mode to consume saved group JSON/PNGs, without rerunning
  scientific analyses. Validate both variants before replacing HTML.
- Before a new deployment, run `meg/tests/test_drsa_group_reports.py` and verify
  the renderer is group-only. Do not reuse an older frozen runtime for a new
  run without checking this policy. Preserve frozen numerical deployments.
- This policy was approved in chat `01a1f428-bd27-73d3-aeb6-6bfcf43122db`
  (19 November 2026), implemented as `c2a3a60` in dad-cluster-analysis, and
  restored in the MNN branch as `fc4aabf`. Deploy a reviewed revision that
  contains this policy and its regression tests; verify the actual checkout
  before submitting jobs. Local worktree paths are not deployment guarantees.
