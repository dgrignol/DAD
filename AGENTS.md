# AGENTS.md

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
