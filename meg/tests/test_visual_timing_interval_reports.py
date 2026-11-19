"""Verify added Stage 3 interval views against independent sample expectations.

Usage from DAD: PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg python3 -B -m pytest \
    -p no:cacheprovider meg/tests/test_visual_timing_interval_reports.py -q

Only synthetic measured events are used. Native 500-Hz samples, nonzero
first_samp and 120-Hz motion frames define known intervals. Tests reject broken
joins and arithmetic, check full-range point retention including sparse/tiny
variance and tails, and verify categorical positions and original bar retention.
No participant signals, shared derivatives or scheduler are touched.
"""

from __future__ import annotations

import copy
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code/preprocessing_basic"))
import visual_timing_report as report


# Explicit sample/frame coordinates provide an independent arithmetic oracle.
def evidence():
    """Make three endpoints sharing one trial, with residuals 2, 2 and 0 ms."""
    rows = []
    for code,digital,edge,frame in ((81,1000,1009,1),(111,1050,1060,13),(112,1100,1110,25)):
        rows.append(dict(event_key=f"sub-11:one.fif:{digital}:{code}",source_fif="one.fif",runtime="v25",
            measurement_status="measured",trigger_code=code,digital_sample_raw=digital,
            diode_edge_sample_raw=edge,first_samp=500,sfreq_hz=500,behavioral_frame=frame,
            lag_ms=(edge-digital)*2,behavioral_block="3",behavioral_run="1",behavioral_trial="7",
            is_replay_trial="False",reported_source_csv="saved.csv",saved_frame_missed_count="0"))
    intervals = [dict(left_event_key=rows[a]["event_key"],right_event_key=rows[b]["event_key"],
                      expected_frame_interval_ms=expected,observed_diode_interval_ms=observed,interval_error_ms=observed-expected)
                 for a,b,expected,observed in ((0,1,100,102),(0,2,200,202),(1,2,100,100))]
    return intervals,rows,{"3":{"refrate":120}}


def test_native_sample_frame_context_and_descriptive_counts():
    inputs = evidence();before = copy.deepcopy(inputs)
    pairs = report.interval_observations(*inputs)
    assert [r["residual_ms"] for r in pairs]==[2,2,0]
    assert [r["recording_time_sec"] for r in pairs]==[1,1,1.1]
    summary = report.interval_statistics(pairs)[0]
    assert summary["n"]==3 and summary["n_trial_contexts"]==1
    assert summary["mean"]==pytest.approx(4/3) and summary["sd"]==pytest.approx(np.sqrt(4/3))
    assert len(report.interval_statistics(pairs,by_transition=True))==3
    assert inputs==before


@pytest.mark.parametrize('problem', ['context','status','rate','sample','frame','duplicate','missing'])
def test_contradictory_interval_evidence_is_rejected(problem):
    intervals,rows,contexts = evidence()
    if problem=='context':rows[1]['behavioral_trial']='8'
    elif problem=='status':rows[1]['measurement_status']='flagged'
    elif problem=='rate':rows[1]['sfreq_hz']=1000
    elif problem=='sample':rows[1]['diode_edge_sample_raw']+=1
    elif problem=='frame':contexts['3']['refrate']=60
    elif problem=='duplicate':intervals.append(copy.deepcopy(intervals[0]))
    elif problem=='missing':rows.pop()
    with pytest.raises(ValueError):report.interval_observations(intervals,rows,contexts)


# Axis inspection verifies the requested encoding without brittle image hashes.
@pytest.mark.parametrize('values', [[0],[0,0,0,0],[.1,.10000000000001,.1],[0,1,2,3],[0,0,0,301]])
def test_violin_retains_every_value_including_constant_and_tail(values,monkeypatch):
    captured=[]
    monkeypatch.setattr(report,'png',lambda fig:captured.append(fig) or b'png')
    rows=[dict(left_event_key=f'sub-11:one.fif:{i}:81',right_event_key=f'sub-11:one.fif:{i+10}:111',
               expected_frame_interval_ms=100,observed_diode_interval_ms=100+v) for i,v in enumerate(values)]
    before=copy.deepcopy(rows)
    report.interval_violin_figure(rows)
    axis=captured[0].axes[1]
    cloud=next(c for c in axis.collections if c.get_label()=='One measured interval')
    offsets=np.asarray(cloud.get_offsets())
    assert len(offsets)==len(values) and offsets[:,1]==pytest.approx(values)
    assert all(-.14<=x<=.14 for x in offsets[:,0])
    assert axis.get_ylim()[0]<=min(values) and axis.get_ylim()[1]>=max(values)
    assert not axis.patches and rows==before
    report.pyplot().close(captured[0])


def test_original_bars_retained_and_new_left_coordinates_match(monkeypatch):
    captured=[]
    monkeypatch.setattr(report,'png',lambda fig:captured.append(fig) or b'png')
    intervals,_,_=evidence()
    report.interval_figure(intervals)
    report.interval_violin_figure(intervals)
    assert len(captured[0].axes[1].patches)==1
    assert captured[0].axes[1].patches[0].get_height()==pytest.approx(4/3)
    old,new=captured[0].axes[0],captured[1].axes[0]
    assert old.get_legend_handles_labels()[1]==new.get_legend_handles_labels()[1]
    np.testing.assert_array_equal(old.collections[0].get_offsets(),new.collections[0].get_offsets())
    for fig in captured:report.pyplot().close(fig)


def test_supplemental_values_category_order_and_stable_jitter(monkeypatch):
    captured=[]
    monkeypatch.setattr(report,'png',lambda fig:captured.append(fig) or b'png')
    pairs=report.interval_observations(*evidence())
    before=copy.deepcopy(pairs)
    files=report.interval_residual_figures(pairs)
    assert set(files)=={'interval_residual_vs_duration.png','interval_residual_vs_time.png','interval_residual_by_pair.png'}
    for figure,expected_x in zip(captured[:2],([100,200,100],[1,1,1.1])):
        points=np.concatenate([np.asarray(c.get_offsets()) for c in figure.axes[0].collections])
        np.testing.assert_allclose(points[:,0],expected_x,rtol=0,atol=1e-12)
        np.testing.assert_allclose(points[:,1],[2,2,0],rtol=0,atol=1e-12)
    pair_axis=captured[2].axes[0]
    assert [t.get_text() for t in pair_axis.get_xticklabels()]==['81->111','81->112','111->112']
    assert list(pair_axis.get_xticks())==[0,1,2]
    assert report.interval_jitter(pairs[0])==report.interval_jitter(copy.deepcopy(pairs[0]))
    assert pairs==before
    for figure in captured:report.pyplot().close(figure)


def test_unavailable_pairs_keep_explanations_without_observations():
    assert report.interval_observations([],[],{})==[]
    files=report.interval_residual_figures([])
    assert len(files)==3 and all(v.startswith(b'\x89PNG') for v in files.values())

