#!/usr/bin/env python3
"""Record reconciled Stage 4 decisions in the existing subject tuning YAML.

Usage after independent inspection and reconciliation
---------------------------------------------------
    python meg/code/preprocessing_basic/review_cleaning.py \
      --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
      --tuning-config /path/to/03_tuning_choices.yaml --stage 4A \
      --reviewer 'User and Codex' --reason 'Reviewed per-run channels and SSS QC'

For 4C add --exclude 0 2 (or --exclude with no numbers for reviewed empty).
4B with no additional high-pass uses no --filter-edge-margin-sec argument.
Spectrum-fit boundary/short-stretch review remains necessary for every source,
including P38 run 1. Only an applied optional high-pass requires an explicit
margin covering its FIR half-support; exposure never excludes epochs. Existing
historical margins/receipts are not rewritten. No processing or submission occurs.

The command binds exact existing report evidence and the ICA model when relevant.
It does not create scientific acceptance merely by being run: use only after
the humans have reconciled the evidence. Formal acceptance remains in the shared
PREPROCESSING_DECISIONS.md. Final completion requires separate applied-output QC.
Under the explicit 2026-11-07 current-annotation policy, the separate
stage4b_adoption.py cohort reconciliation records processing acceptance only
after validating all operative outputs. Those runs retain unreviewed manual
status; a bound adoption plan supplies independent processing provenance.
"""

import argparse
import json
import math
from pathlib import Path

from continuous_cleaning import file_hash, require_review, stage_configuration
from pipeline_config import load_and_resolve_config
from preprocessing_stage_utils import update_tuning_file, utc_now


def main(argv=None):
    """Bind an explicit human decision; no implicit defaults represent review."""
    p = argparse.ArgumentParser(description=__doc__)
    for flag in ("config", "subject-config", "tuning-config", "reviewer", "reason"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--site-config")
    p.add_argument("--stage", choices=("4A", "4B", "4C"), required=True)
    p.add_argument("--exclude", nargs="*", type=int)
    p.add_argument("--filter-edge-margin-sec", type=float)
    args = p.parse_args(argv)
    cfg, *_ = load_and_resolve_config(args.config, args.subject_config, args.tuning_config if Path(args.tuning_config).exists() else None,
                                      site_config_path=args.site_config, required_sections=())
    report_dir = Path(cfg["resolved_paths"]["derivatives_reports_dir"])
    name = {"4A": "03_A_MaxFilt", "4B": "03_B_ArtefactAnnotation", "4C": "03_C_ICA"}[args.stage]
    report = json.loads((report_dir / f"{name}_review_report.json").read_text())
    if report.get("subject_bids_id") != cfg["layout"]["subject_bids_id"] or not report.get("evidence_sha256"):
        raise ValueError("Report lacks current Stage 4 evidence for this participant.")
    if report.get("configuration_sha256") != stage_configuration(cfg, name, include_edge=False):
        raise ValueError("Configuration differs from the report; regenerate affected review evidence first.")
    # Validate accepted files before creating a receipt, using the same checker
    # that downstream stages run. No acceptance is inferred from a filename.
    trial = dict(cfg)
    trial["stage4_review"] = {**cfg.get("stage4_review", {}), name: {
        "status": "accepted", "evidence_sha256": report["evidence_sha256"],
        "configuration_sha256": stage_configuration(cfg, name)}}
    require_review(trial, name, report_dir)
    if args.stage != "4A":
        require_review(cfg, "03_A_MaxFilt", report_dir)
    if args.stage == "4C":
        require_review(cfg, "03_B_ArtefactAnnotation", report_dir)
        if report.get("component_review_coverage", {}).get("complete") is False:
            p.error("Final-model component diagnostics are incomplete; regenerate them before review acceptance.")
        if args.exclude is None:
            p.error("4C requires --exclude; an explicitly empty list is valid.")
        if any(i < 0 or i >= report["n_components_fitted"] for i in args.exclude):
            p.error("Component outside fitted range.")
        model = Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]) / "ica_models" / Path(report["prefit_model"]).name
        if file_hash(model) != report["prefit_model_sha256"]:
            raise ValueError("Model changed after the component review report.")
        update_tuning_file(path=args.tuning_config, raw_root=cfg["resolved_paths"]["data_root"],
                           subject_id=cfg["subject"]["id"], section="ica",
                           values={"manual_exclude": sorted(set(args.exclude)), "reviewed_model_sha256": file_hash(model)})
    if args.stage == "4B":
        margin = args.filter_edge_margin_sec
        required = max((row["final_filter_preview"]["filter"]["half_support_sec"]
                        for row in report["runs"] if row["final_filter_preview"]["filter"] is not None), default=0.)
        if required and (margin is None or not math.isfinite(margin) or margin < required):
            p.error(f"Explicit informational FIR exposure margin must be at least {required:g} s; review exposure first.")
        if not required and margin is not None:
            p.error("No high-pass is applied: omit --filter-edge-margin-sec; spectrum-fit has no FIR exposure margin.")
        if report.get("annotation_review_policy") == "current_annotations_defer_detailed_review_2026-11-07":
            from stage4b_adoption import require_operative_adoption
            require_operative_adoption(cfg, report, require_publication=True)
        elif not all(row.get("manual_eye_review", {}).get("status") == "reviewed" for row in report["runs"]):
            p.error("Every acquisition needs an explicit ocular review or validated explicit-policy adoption.")
        if any(row.get("edf_blink_candidates", {}).get("review_status") == "review_required" for row in report["runs"]):
            p.error("New/changed EDF guidance needs explicit complete-source ocular review before 4B acceptance.")
        if required:
            update_tuning_file(path=args.tuning_config, raw_root=cfg["resolved_paths"]["data_root"],
                               subject_id=cfg["subject"]["id"], section="preprocessing",
                               values={"filter_edge_margin_sec": margin})
    # Reload deliberate tuning updates before recording the accepted configuration.
    cfg, *_ = load_and_resolve_config(args.config, args.subject_config,
        args.tuning_config if Path(args.tuning_config).exists() else None,
        site_config_path=args.site_config, required_sections=())
    receipt = {"status": "accepted", "configuration_sha256": stage_configuration(cfg, name), "evidence_sha256": report["evidence_sha256"],
               "reviewer": args.reviewer, "reason": args.reason, "reviewed_utc": utc_now()}
    if args.stage == "4B" and report.get("annotation_review_policy"):
        receipt.update(acceptance_kind="explicit_annotation_policy_for_processing", policy=report["annotation_review_policy"])
    update_tuning_file(path=args.tuning_config, raw_root=cfg["resolved_paths"]["data_root"],
                       subject_id=cfg["subject"]["id"], section="stage4_review", values={name: receipt})
    print(f"Recorded {args.stage} decision; formal shared acceptance and applied-output review remain explicit.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
