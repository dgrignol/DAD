#!/usr/bin/env python3
"""
run_pipeline.py

Purpose
Provide one command-line entrypoint for configuration resolution, automatic
dataset preflight, durable source lineage, provenance snapshots, and execution
of pipeline stages that have a stable noninteractive CLI.

Current scope
The stable production stages are trigger extraction, Maxwell filtering,
physiological artifact annotation, manual ICA, epoching, and dRSA neural
export. MATLAB dRSA is registered as an optional downstream analysis stage.
The automatic MEGNet ICA branch is deprecated experimental code and is
intentionally not registered.

Usage examples
1) Preflight only:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/run_pipeline.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --preflight-only

   Add `--write-inventory` when separate inventory JSON/Markdown files are
   wanted. Without it, the preflight result is kept in the pipeline run log.

2) Preview commands without executing trigger extraction:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/run_pipeline.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --stages trigger_extraction --dry-run

3) Run trigger extraction:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/run_pipeline.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --stages trigger_extraction

   The trigger stage writes both the canonical event/alignment derivatives and
   the legacy deep Markdown/CSV/plot report bundle. Add
   `--disable-deep-trigger-report` only when those plots are intentionally not
   required.

4) Run the first-pass production preprocessing review:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/run_pipeline.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --stages maxwell_filter artifact_annotation ica_manual \
     --processing-pass review

5) Apply accepted tuning:
   python .../run_pipeline.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml \
     --tuning-config .../03_tuning_choices.yaml \
     --stages maxwell_filter artifact_annotation ica_manual \
     --processing-pass apply --overwrite

6) Resume interrupted preprocessing stages. Existing derivatives are
   validated before they are skipped:

   python .../run_pipeline.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml \
     --tuning-config .../03_tuning_choices.yaml \
     --stages maxwell_filter artifact_annotation ica_manual \
     --processing-pass review --resume

7) Run the noninteractive review pass on the CIMeC cluster:

   python .../run_pipeline.py --config .../pipeline_default.yaml \
     --site-config .../pipeline_cluster.yaml \
     --subject-config .../sub-02.yaml \
     --stages maxwell_filter artifact_annotation ica_manual \
     --processing-pass review --resume

   The SGE array wrapper and launch guide live in ``meg/cluster``.

8) Run downstream epoching and dRSA export after manual-ICA apply completion:

   python .../run_pipeline.py --config .../pipeline_default.yaml \
     --subject-config .../sub-04.yaml \
     --stages epoching export_drsa

Data flow
configuration -> resolved-config snapshot + source manifest -> automatic
inventory/preflight -> selected stage subprocesses -> pipeline run log.

Safety
The runner validates every infrastructure output against the configured raw
root. Individual stages remain responsible for applying the same guard to
their own outputs.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from time import monotonic
from typing import Any, Mapping

from dataset_inventory import (
    build_inventory_payload,
    build_source_manifest_payload,
    write_inventory_reports,
    write_source_manifest,
)
from pipeline_config import (
    assert_output_path_outside_raw,
    load_and_resolve_config,
    write_resolved_config,
)


SUPPORTED_STAGES = (
    "trigger_extraction",
    "maxwell_filter",
    "artifact_annotation",
    "ica_manual",
    "epoching",
    "export_drsa",
    "matlab_drsa",
)


def _iso_utc_now() -> str:
    """Return a UTC ISO8601 timestamp with explicit Z suffix."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def build_stage_commands(
    *,
    config_path: Path,
    site_config_path: Path | None = None,
    subject_config_path: Path | None,
    source_manifest_path: Path,
    disable_deep_trigger_report: bool,
    tuning_config_path: Path | None = None,
    processing_pass: str = "review",
    interactive: bool = False,
    no_display: bool = False,
    decisions_out: Path | None = None,
    overwrite: bool = False,
    resume: bool = False,
) -> dict[str, list[str]]:
    """
    Build subprocess commands for currently supported pipeline stages.

    Registered stages have stable configuration-driven CLIs. Interactivity is
    opt-in; batch mode remains the default.
    """
    script_dir = Path(__file__).resolve().parent
    trigger_command = [
        sys.executable,
        str(script_dir / "02_trigger_extraction_custom.py"),
        "--config",
        str(config_path),
        "--source-manifest",
        str(source_manifest_path),
        "--print-summary",
    ]
    if subject_config_path is not None:
        trigger_command.extend(["--subject-config", str(subject_config_path)])
    if site_config_path is not None:
        trigger_command.extend(["--site-config", str(site_config_path)])
    if disable_deep_trigger_report:
        trigger_command.append("--disable-deep-trigger-report")

    commands = {"trigger_extraction": trigger_command}
    stage_scripts = {
        "maxwell_filter": "03_A_MaxFilt.py",
        "artifact_annotation": "03_B_ArtefactAnnotation.py",
        "ica_manual": "03_C_ICA.py",
    }
    for stage_name, script_name in stage_scripts.items():
        command = [
            sys.executable,
            str(script_dir / script_name),
            "--config",
            str(config_path),
            "--source-manifest",
            str(source_manifest_path),
            "--processing-pass",
            processing_pass,
            "--print-summary",
        ]
        if subject_config_path is not None:
            command.extend(["--subject-config", str(subject_config_path)])
        if site_config_path is not None:
            command.extend(["--site-config", str(site_config_path)])
        if tuning_config_path is not None:
            command.extend(["--tuning-config", str(tuning_config_path)])
        if interactive:
            command.append("--interactive")
        if no_display:
            command.append("--no-display")
        if decisions_out is not None:
            command.extend(["--decisions-out", str(decisions_out)])
        if overwrite:
            command.append("--overwrite")
        if resume:
            command.append("--resume")
        commands[stage_name] = command

    # Downstream scientific stages consume the completed manifest lineage and
    # do not expose Stage-03 review/apply or interactive decision controls.
    downstream_scripts = {
        "epoching": "04_epoching.py",
        "export_drsa": "06_export_neural_for_drsa.py",
        "matlab_drsa": "07_run_matlab_drsa.py",
    }
    for stage_name, script_name in downstream_scripts.items():
        command = [
            sys.executable,
            str(script_dir / script_name),
            "--config",
            str(config_path),
            "--source-manifest",
            str(source_manifest_path),
            "--print-summary",
        ]
        if subject_config_path is not None:
            command.extend(["--subject-config", str(subject_config_path)])
        if site_config_path is not None:
            command.extend(["--site-config", str(site_config_path)])
        if tuning_config_path is not None:
            command.extend(["--tuning-config", str(tuning_config_path)])
        if overwrite:
            command.append("--overwrite")
        if resume:
            command.append("--resume")
        commands[stage_name] = command
    return commands


def _write_run_log(
    path: str | Path,
    payload: Mapping[str, Any],
    raw_root: str | Path,
) -> Path:
    """Write the machine-readable orchestration log outside the raw tree."""
    resolved = assert_output_path_outside_raw(path, raw_root)
    resolved.parent.mkdir(parents=True, exist_ok=True)
    resolved.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return resolved


def preserve_durable_manifest_sections(
    new_manifest: dict[str, Any],
    existing_manifest_path: Path,
) -> dict[str, Any]:
    """
    Preserve evidence-derived manifest sections across orchestration runs.

    Source discovery is refreshed from configuration, while trigger-derived
    associations and scientific derivative lineage remain durable.
    """
    if not existing_manifest_path.is_file():
        return new_manifest
    try:
        existing = json.loads(
            existing_manifest_path.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError):
        return new_manifest
    for key in ("trigger_derived_associations", "derivative_lineage"):
        if key in existing:
            new_manifest[key] = existing[key]
    return new_manifest


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the orchestration CLI."""
    parser = argparse.ArgumentParser(description="Run the DAD MEG pipeline.")
    parser.add_argument("--config", required=True)
    parser.add_argument(
        "--site-config",
        default=None,
        help="Optional machine/site path override YAML.",
    )
    parser.add_argument("--subject-config", default=None)
    parser.add_argument(
        "--stages",
        nargs="+",
        choices=SUPPORTED_STAGES,
        default=["trigger_extraction"],
    )
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--skip-preflight", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument("--disable-deep-trigger-report", action="store_true")
    parser.add_argument(
        "--processing-pass",
        choices=("review", "apply"),
        default="review",
        help="Review generates QC/tuning evidence; apply requires accepted ICA choices.",
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Prompt for bad-channel, muscle-threshold, and ICA decisions.",
    )
    parser.add_argument(
        "--no-display",
        action="store_true",
        help="With --interactive, save plots without opening GUI windows.",
    )
    parser.add_argument(
        "--tuning-config",
        default=None,
        help="Accepted review decisions YAML used for an apply pass.",
    )
    parser.add_argument(
        "--decisions-out",
        default=None,
        help="Interactive decisions YAML; defaults to 03_tuning_choices.yaml.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Intentionally replace existing scientific stage derivatives.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Resume supported preprocessing stages by validating and skipping "
            "existing derivatives. Mutually exclusive with --overwrite."
        ),
    )
    parser.add_argument("--resolved-config-out", default=None)
    parser.add_argument(
        "--source-manifest",
        default=None,
        help=(
            "Durable source-lineage JSON. Defaults to source_manifest.json in "
            "the subject derivatives report directory."
        ),
    )
    parser.add_argument(
        "--write-inventory",
        action="store_true",
        help=(
            "Write default dataset_inventory.json/.md files. Preflight still "
            "runs without this flag and is summarized in the pipeline run log."
        ),
    )
    parser.add_argument("--inventory-json", default=None)
    parser.add_argument("--inventory-md", default=None)
    parser.add_argument("--run-log", default=None)
    return parser


def main() -> int:
    """Resolve inputs, run preflight, execute selected stages, and log results."""
    args = build_arg_parser().parse_args()
    if args.resume and args.overwrite:
        print(
            "Pipeline configuration error: --resume and --overwrite are "
            "mutually exclusive.",
            file=sys.stderr,
        )
        return 2
    try:
        required_sections = []
        if "trigger_extraction" in args.stages:
            required_sections.append("trigger_parser")
        if "maxwell_filter" in args.stages:
            required_sections.append("maxwell_filter")
        if "artifact_annotation" in args.stages:
            required_sections.append("artifact_annotation")
        if "ica_manual" in args.stages:
            required_sections.append("ica")
        if (
            "epoching" in args.stages
            or "export_drsa" in args.stages
            or "matlab_drsa" in args.stages
        ):
            required_sections.append("epoching")
        if "export_drsa" in args.stages or "matlab_drsa" in args.stages:
            required_sections.append("export")
        if "matlab_drsa" in args.stages:
            required_sections.append("matlab_drsa")
        cfg, warnings, config_path, subject_config_path = load_and_resolve_config(
            args.config,
            args.subject_config,
            args.tuning_config,
            site_config_path=args.site_config,
            required_sections=tuple(required_sections),
        )
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Pipeline configuration error: {exc}", file=sys.stderr)
        return 2

    resolved = cfg["resolved_paths"]
    raw_root = resolved["data_root"]
    reports_dir = Path(resolved["derivatives_reports_dir"])
    cleaning_only = bool(args.stages) and set(args.stages) <= {"maxwell_filter", "artifact_annotation", "ica_manual"}
    snapshot_prefix = "stage4_" if cleaning_only else ""
    resolved_config_out = (
        Path(args.resolved_config_out).expanduser().resolve()
        if args.resolved_config_out
        else reports_dir / f"{snapshot_prefix}resolved_config.json"
    )
    inventory_json = None
    if args.inventory_json:
        inventory_json = Path(args.inventory_json).expanduser().resolve()
    elif args.write_inventory:
        inventory_json = reports_dir / "dataset_inventory.json"

    inventory_md = None
    if args.inventory_md:
        inventory_md = Path(args.inventory_md).expanduser().resolve()
    elif args.write_inventory:
        inventory_md = reports_dir / "dataset_inventory.md"
    run_log_path = (
        Path(args.run_log).expanduser().resolve()
        if args.run_log
        else reports_dir / f"{snapshot_prefix}pipeline_run.json"
    )
    source_manifest_path = (
        Path(args.source_manifest).expanduser().resolve()
        if args.source_manifest
        else reports_dir / "source_manifest.json"
    )
    tuning_config_path = (
        Path(args.tuning_config).expanduser().resolve()
        if args.tuning_config
        else None
    )
    site_config_path = (
        Path(args.site_config).expanduser().resolve()
        if args.site_config
        else None
    )
    decisions_out = (
        Path(args.decisions_out).expanduser().resolve()
        if args.decisions_out
        else reports_dir / "03_tuning_choices.yaml"
    )

    run_payload: dict[str, Any] = {
        "started_utc": _iso_utc_now(),
        "subject_bids_id": cfg["layout"]["subject_bids_id"],
        "session_bids_id": cfg["layout"]["session_bids_id"],
        "dry_run": bool(args.dry_run),
        "processing_pass": args.processing_pass,
        "interactive": bool(args.interactive),
        "tuning_config": str(tuning_config_path) if tuning_config_path else None,
        "site_config": str(site_config_path) if site_config_path else None,
        "decisions_out": str(decisions_out),
        "overwrite": bool(args.overwrite),
        "resume": bool(args.resume),
        "preflight_skipped": bool(args.skip_preflight),
        "selected_stages": list(args.stages),
        "warnings": warnings,
        "stages": [],
    }

    try:
        written_config = write_resolved_config(resolved_config_out, cfg)
        run_payload["resolved_config"] = str(written_config)

        # The source manifest is durable provenance. It is independent of the
        # optional inventory reports and never asserts positional alignment
        # between MEG acquisition runs and behavioral blocks.
        if cleaning_only:
            if not source_manifest_path.is_file():
                raise RuntimeError("Stage 4 requires the existing accepted source manifest.")
            run_payload["source_manifest"] = str(source_manifest_path)
        else:
            source_manifest = build_source_manifest_payload(cfg, config_path, subject_config_path, warnings)
            source_manifest = preserve_durable_manifest_sections(source_manifest, source_manifest_path)
            written_manifest = write_source_manifest(source_manifest, source_manifest_path)
            run_payload["source_manifest"] = str(written_manifest)

        if not args.skip_preflight and not cleaning_only:
            inventory = build_inventory_payload(
                cfg,
                config_path,
                subject_config_path,
                warnings,
            )
            inventory_outputs = write_inventory_reports(
                inventory,
                out_json=inventory_json,
                out_md=inventory_md,
            )
            run_payload["inventory_status"] = inventory["status"]
            run_payload["inventory_outputs"] = inventory_outputs
            run_payload["inventory_summary"] = {
                "counts": inventory["counts"],
                "checks": inventory["checks"],
                "warnings": inventory["validation"]["warnings"],
                "errors": inventory["validation"]["errors"],
            }
            if not inventory["validation"]["is_valid"]:
                raise RuntimeError(
                    "Dataset preflight failed: "
                    + "; ".join(inventory["validation"]["errors"])
                )
        else:
            run_payload["inventory_status"] = "reused_accepted_stage1" if cleaning_only else "skipped"

        commands = build_stage_commands(
            config_path=config_path,
            site_config_path=site_config_path,
            subject_config_path=subject_config_path,
            source_manifest_path=source_manifest_path,
            disable_deep_trigger_report=args.disable_deep_trigger_report,
            tuning_config_path=tuning_config_path,
            processing_pass=args.processing_pass,
            interactive=args.interactive,
            no_display=args.no_display,
            decisions_out=decisions_out,
            overwrite=args.overwrite,
            resume=args.resume,
        )
        for stage in args.stages:
            command = commands[stage]
            stage_record: dict[str, Any] = {
                "stage": stage,
                "command": command,
                "command_display": shlex.join(command),
                "started_utc": _iso_utc_now(),
            }
            run_payload["stages"].append(stage_record)
            print(f"[{stage}] {stage_record['command_display']}")

            if args.preflight_only or args.dry_run:
                stage_record["status"] = "not_run"
                stage_record["reason"] = (
                    "preflight_only" if args.preflight_only else "dry_run"
                )
                continue

            start = monotonic()
            completed = subprocess.run(command, check=False)
            stage_record["duration_sec"] = monotonic() - start
            stage_record["return_code"] = int(completed.returncode)
            stage_record["status"] = (
                "completed" if completed.returncode == 0 else "failed"
            )
            if completed.returncode != 0 and not args.continue_on_error:
                break
    except Exception as exc:  # pylint: disable=broad-except
        run_payload["orchestration_error"] = str(exc)
        print(f"Pipeline error: {exc}", file=sys.stderr)

    run_payload["finished_utc"] = _iso_utc_now()
    failed = bool(run_payload.get("orchestration_error")) or any(
        item.get("status") == "failed" for item in run_payload["stages"]
    )
    run_payload["status"] = "failed" if failed else "completed"
    try:
        written_log = _write_run_log(run_log_path, run_payload, raw_root)
        print(f"Pipeline run log: {written_log}")
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Could not write pipeline run log: {exc}", file=sys.stderr)
        return 2
    return 2 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
