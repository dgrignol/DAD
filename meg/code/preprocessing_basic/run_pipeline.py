#!/usr/bin/env python3
"""
run_pipeline.py

Purpose
Provide one command-line entrypoint for configuration resolution, automatic
dataset preflight, provenance snapshots, and execution of pipeline stages that
have a stable noninteractive CLI.

Current scope
Only `trigger_extraction` is registered. The remaining preprocessing scripts
are intentionally not exposed until their hardcoded paths and subject settings
are refactored. Future stages should be added to `build_stage_commands()` after
they accept the same config/subject-config contract.

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

Data flow
configuration -> resolved-config snapshot -> automatic inventory/preflight ->
selected stage subprocesses -> pipeline run log.

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

from dataset_inventory import build_inventory_payload, write_inventory_reports
from pipeline_config import (
    assert_output_path_outside_raw,
    load_and_resolve_config,
    write_resolved_config,
)


SUPPORTED_STAGES = ("trigger_extraction",)


def _iso_utc_now() -> str:
    """Return a UTC ISO8601 timestamp with explicit Z suffix."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def build_stage_commands(
    *,
    config_path: Path,
    subject_config_path: Path | None,
    disable_deep_trigger_report: bool,
) -> dict[str, list[str]]:
    """
    Build subprocess commands for currently supported pipeline stages.

    New stages belong here only after they are configuration-driven and
    noninteractive. This keeps `--dry-run` honest about what can be automated.
    """
    script_dir = Path(__file__).resolve().parent
    trigger_command = [
        sys.executable,
        str(script_dir / "02_trigger_extraction_custom.py"),
        "--config",
        str(config_path),
        "--print-summary",
    ]
    if subject_config_path is not None:
        trigger_command.extend(["--subject-config", str(subject_config_path)])
    if disable_deep_trigger_report:
        trigger_command.append("--disable-deep-trigger-report")
    return {"trigger_extraction": trigger_command}


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


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the orchestration CLI."""
    parser = argparse.ArgumentParser(description="Run the DAD MEG pipeline.")
    parser.add_argument("--config", required=True)
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
    parser.add_argument("--resolved-config-out", default=None)
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
    try:
        cfg, warnings, config_path, subject_config_path = load_and_resolve_config(
            args.config,
            args.subject_config,
            required_sections=("trigger_parser",),
        )
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Pipeline configuration error: {exc}", file=sys.stderr)
        return 2

    resolved = cfg["resolved_paths"]
    raw_root = resolved["data_root"]
    reports_dir = Path(resolved["derivatives_reports_dir"])
    resolved_config_out = (
        Path(args.resolved_config_out).expanduser().resolve()
        if args.resolved_config_out
        else reports_dir / "resolved_config.json"
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
        else reports_dir / "pipeline_run.json"
    )

    run_payload: dict[str, Any] = {
        "started_utc": _iso_utc_now(),
        "subject_bids_id": cfg["layout"]["subject_bids_id"],
        "session_bids_id": cfg["layout"]["session_bids_id"],
        "dry_run": bool(args.dry_run),
        "preflight_skipped": bool(args.skip_preflight),
        "selected_stages": list(args.stages),
        "warnings": warnings,
        "stages": [],
    }

    try:
        written_config = write_resolved_config(resolved_config_out, cfg)
        run_payload["resolved_config"] = str(written_config)

        if not args.skip_preflight:
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
            run_payload["inventory_status"] = "skipped"

        commands = build_stage_commands(
            config_path=config_path,
            subject_config_path=subject_config_path,
            disable_deep_trigger_report=args.disable_deep_trigger_report,
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
