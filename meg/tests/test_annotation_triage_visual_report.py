"""Validate portable retention books with synthetic saved trial/count snapshots.

Usage:
    python -B -m pytest -p no:cacheprovider meg/tests/test_annotation_triage_visual_report.py -q

Exercise real plot rendering, embedded HTML images, reconciliation with trial
ledgers, all 120 identity cells, zero denominators, incomplete-input handling,
annotation-only baseline comparisons, immutable publication and source protection.
No actual participant files, MEG recordings, detectors or scheduler are used.
"""

import csv
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
import annotation_triage as triage
import render_annotation_triage_report as visual


# -----------------------------------------------------------------------------
# Small trial-level evidence includes overlaps, absent run conditions and errors
# -----------------------------------------------------------------------------

def snapshot(tmp_path, label="before", *, revised=False, zero=False, incomplete=False):
    """Write the same table contract as the counter without reading any FIF."""
    folder = tmp_path / "annotation_triage/reports" / label
    folder.mkdir(parents=True)
    cfg = {"subject_profile": {"identity_range": [1,40]},
        "epoching": {"event_id_map": {"visible":71,"occluded":81,"deviant":91}}}
    trials=[]
    for code in (71,81,91):
        for identity in range(1,41):
            bad = [name for divisor,name in ((3,"BAD_muscle"),(5,"BAD_blink_auto")) if identity % divisor==0]
            if revised and identity in (3,6): bad=[]
            if revised and identity==7: bad=["BAD_blink_auto"]
            if zero: bad=["BAD_muscle"]
            boundary=identity==40
            kept=not bad and not boundary
            trials.append({"subject_bids_id":"sub-06","source_fif":"run01.fif","trial_key":f"c{code}-i{identity}",
                "condition_code":code,"trial_sequence_identity":identity,"analysis_eligible":True,
                "annotation_survives":kept,"incomplete_window":boundary,"bad_categories":";".join(bad),
                "unscored_groups":"eye_movement","epoch_start_raw":identity*1000,"epoch_stop_raw_exclusive":identity*1000+318,
                "triage_outcome":"incomplete_window" if boundary else "overlaps_BAD" if bad else "annotation_survivor"})
    # A partial acquisition can truly have zero trials for some conditions.
    trials.append({**trials[0],"source_fif":"run02.fif","trial_key":"partial","analysis_eligible":False,
        "annotation_survives":False,"trial_sequence_identity":None,"incomplete_window":False,
        "triage_outcome":"excluded_by_accepted_trial_policy","epoch_start_raw":"","epoch_stop_raw_exclusive":""})
    conditions,identities,runs=triage.summarize(trials,cfg,["run01.fif","run02.fif"])
    for table in (conditions,identities,runs):
        for row in table: row.update(subject_bids_id="sub-06",input_status="counted")
    summary={"subject_bids_id":"sub-06","status":"counted",
        "minimum_condition_survivors":min(r["annotation_survivors"] for r in conditions),
        "total_annotation_survivors":sum(r["annotation_survivors"] for r in conditions),
        "missing_identity_cells":sum(r["annotation_survivors"]==0 for r in identities),
        "runs_missing_selections":0,"survivors_with_unscored_ocular":sum(t["annotation_survives"] for t in trials),
        "masks":[{"method":"provisional_automatic"}]}
    payload={"schema":"annotation-triage-v1","purpose":"preliminary_annotation_retention_only",
        "status":"completed","subjects":[summary],"errors":[]}
    if incomplete:
        payload["status"]="incomplete_inputs"
        payload["subjects"].append({"subject_bids_id":"sub-07","status":"incomplete_inputs"})
        payload["errors"].append({"subject_bids_id":"sub-07","source_fif":"bad.fif","error":"stale <model>"})
    (folder/"report.json").write_text(json.dumps(payload))
    for name,rows in zip(visual.TABLES,(trials,conditions,identities,runs)):
        (folder/f"{name}.csv").write_text(triage.csv_content(rows))
    return folder


def change_table(path, mutate):
    """Introduce an inconsistency while retaining well-formed CSV syntax."""
    rows=visual.read_table(path); mutate(rows)
    path.write_text(triage.csv_content(rows))


# -----------------------------------------------------------------------------
# Numerical meanings and annotation-only comparison contracts
# -----------------------------------------------------------------------------

def test_all_identity_bars_reconcile_and_share_scale_with_true_zeros(tmp_path):
    loaded=visual.load_snapshot(snapshot(tmp_path,zero=True))
    part=loaded["participants"]["sub-06"]
    assert len(part["identity_retention"])==120
    fig=visual.identity_figure(part)
    assert len(fig.axes)==3
    assert len({ax.get_ylim() for ax in fig.axes})==1
    for ax in fig.axes:
        assert len(ax.containers[0])==40
        assert sum(p.get_height() for p in ax.containers[0])==0
        assert sum(p.get_height() for p in ax.containers[1])==40
    assert loaded["ranked"]==["sub-06"]
    rows=visual.summary_rows(loaded)
    assert all(r["empty_identities"]==40 and r["survivors"]==0 for r in rows)
    assert np.isnan(visual.percent(0,0))


@pytest.mark.parametrize("table,mutate",[
    ("identity_retention",lambda r:r.pop()),
    ("identity_retention",lambda r:r[0].update(annotation_survivors=999)),
    ("condition_retention",lambda r:r[0].update(eligible=999)),
    ("trial_retention",lambda r:r[0].update(analysis_eligible="")),
    ("trial_retention",lambda r:r[0].update(annotation_survives="False")),
    ("trial_retention",lambda r:r[0].update(triage_outcome="overlaps_BAD")),
])
def test_malformed_or_inconsistent_cells_fail_before_plotting(tmp_path,table,mutate):
    folder=snapshot(tmp_path)
    change_table(folder/f"{table}.csv",mutate)
    with pytest.raises(ValueError): visual.load_snapshot(folder)


def test_overlap_reasons_are_not_a_unique_loss_sum(tmp_path):
    loaded=visual.load_snapshot(snapshot(tmp_path))
    trials=loaded["participants"]["sub-06"]["trial_retention"]
    counts=visual.reason_counts(trials)
    unique=sum(t["analysis_eligible"] and not t["annotation_survives"] for t in trials)
    assert sum(counts.values())>unique
    assert counts["incomplete_window"]==3


def test_comparison_accepts_annotation_changes_but_rejects_changed_windows(tmp_path):
    before=visual.load_snapshot(snapshot(tmp_path))
    after_dir=snapshot(tmp_path,"after",revised=True)
    after=visual.load_snapshot(after_dir)
    assert visual.compare_snapshots(after,before)=={"sub-06"}
    assert after["summaries"]["sub-06"]["total_annotation_survivors"]>before["summaries"]["sub-06"]["total_annotation_survivors"]
    change_table(after_dir/"trial_retention.csv",lambda rows:rows[0].update(epoch_start_raw="999999"))
    with pytest.raises(ValueError,match="windows changed"):
        visual.compare_snapshots(visual.load_snapshot(after_dir),before)


# -----------------------------------------------------------------------------
# Real self-contained HTML publication and preservation of inputs/prior reports
# -----------------------------------------------------------------------------

def test_render_embeds_all_figures_and_comparison_preserves_sources(tmp_path):
    before=snapshot(tmp_path,incomplete=True)
    after=snapshot(tmp_path,"after",revised=True,incomplete=True)
    original={p:visual.sha(p) for folder in (before,after) for p in folder.iterdir()}
    output=visual.render_report(after,label="books",compare_to=before)
    books=list(output.glob("*.html"))
    assert {p.name for p in books}=={"index.html","sub-06.html","sub-07.html"}
    for path in books:
        text=path.read_text()
        assert '<html lang="en">' in text
        assert 'src="figures/' not in text and '<script src=' not in text
        assert 'http://' not in text and 'https://' not in text
    assert (output/"index.html").read_text().count('src="data:image/png;base64,')==5
    assert (output/"sub-06.html").read_text().count('src="data:image/png;base64,')==3
    assert "Hollow points" in (output/"sub-06.html").read_text()
    error_book=(output/"sub-07.html").read_text()
    assert "unranked" in error_book and "stale &lt;model&gt;" in error_book
    changes=visual.read_table(output/"changes.csv")
    assert len(changes)==3 and all(int(r["change"])>0 for r in changes)
    assert all(visual.sha(p)==expected for p,expected in original.items())
    manifest=json.loads((output/"manifest.json").read_text())
    assert manifest["compared_participants"]==["sub-06"]
    assert all(visual.sha(output/p)==expected for p,expected in manifest["output_sha256"].items())
    with pytest.raises(FileExistsError): visual.render_report(after,label="books")


def test_only_input_errors_still_get_an_honest_report(tmp_path):
    folder=snapshot(tmp_path,incomplete=True)
    payload=json.loads((folder/"report.json").read_text())
    payload["subjects"]=payload["subjects"][1:]
    (folder/"report.json").write_text(json.dumps(payload))
    for name in visual.TABLES: (folder/f"{name}.csv").write_text("")
    output=visual.render_report(folder,label="errors")
    assert "counted: 0/1" in (output/"index.html").read_text()
    assert not list((output/"figures").iterdir())


def test_changed_source_during_render_never_publishes(tmp_path,monkeypatch):
    folder=snapshot(tmp_path)
    original=visual.cohort_book
    def changing(*args,**kwargs):
        original(*args,**kwargs)
        (folder/"report.json").write_text((folder/"report.json").read_text()+"\n")
    monkeypatch.setattr(visual,"cohort_book",changing)
    with pytest.raises(RuntimeError,match="changed during rendering"):
        visual.render_report(folder,label="stale")
    assert not (folder.parent.parent/"visual_reports/stale").exists()
    assert not list((folder.parent.parent/"visual_reports").glob(".building-*"))


def test_parallel_participant_books_are_identical_to_serial_output(tmp_path):
    """Compare real two-participant process rendering with the serial reference."""
    folder=snapshot(tmp_path)
    report=json.loads((folder/"report.json").read_text())
    other=json.loads(json.dumps(report["subjects"][0]));other["subject_bids_id"]="sub-07"
    report["subjects"].append(other)
    (folder/"report.json").write_text(json.dumps(report))
    for name in visual.TABLES:
        rows=visual.read_table(folder/f"{name}.csv")
        rows += [{**r,"subject_bids_id":"sub-07"} for r in rows]
        (folder/f"{name}.csv").write_text(triage.csv_content(rows))
    serial=visual.render_report(folder,label="serial",jobs=1)
    concurrent=visual.render_report(folder,label="parallel",jobs=2)
    for path in serial.rglob("*"):
        if path.is_file() and path.name!="manifest.json":
            assert visual.sha(path)==visual.sha(concurrent/path.relative_to(serial))
