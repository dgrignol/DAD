"""Regression tests of the native short-fragment line-removal guard.

Usage: python -m pytest meg/tests/test_epoch_short_fragments.py -q
Only synthetic in-memory data are processed; participant inputs stay read-only.
Tests exercise real spectrum fitting, empty/single-sample segments, BAD and EDGE
markers, nonzero sample origins, annotation preservation and actual epoch parity.
The pinned cluster stack is MNE 1.10.2, NumPy 1.26.4 and SciPy 1.13.1.
"""
from datetime import datetime,timezone
from pathlib import Path

import mne
import numpy as np
import pytest
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))

SETTINGS={'freqs':[50.,100.],'method':'spectrum_fit','filter_length':'10s'}


def guarded_filter(raw):
    """Exercise the actual production guard and its absolute-sample QC records."""
    from epoch_input_preparation import filter_continuous_for_epoching
    settings={'notch_freqs_hz':[50.,100.], 'notch_method':'spectrum_fit',
        'spectrum_fit_window':'10s', 'highpass_hz':None, 'lowpass_hz':None,
        'resample_sfreq_hz':None}
    summary=filter_continuous_for_epoching(raw,settings,source_filename='synthetic_run.fif')['spectrum_fit']
    assert summary['minimum_segment_samples']==9
    skipped=summary['skipped_fragments']
    assert summary['skipped_seconds']==sum(row['n_samples'] for row in skipped)/raw.info['sfreq']
    return [(row['first_sample_raw']-raw.first_samp,
             row['stop_sample_raw_exclusive']-raw.first_samp) for row in skipped]


@pytest.mark.parametrize('length',list(range(17)))
@pytest.mark.parametrize('date',[False,True])
def test_guard_matches_raw_reference_and_epochs(length,date):
    """Compare with independent Raw-level exclusion of the tiny filter island."""
    sfreq=1000.; size=16000
    rng=np.random.default_rng(19)
    values=rng.normal(size=(3,size));values[0]*=1e-13;values[2]=np.arange(size)%250
    times=np.arange(size)/sfreq
    values[0]+=2e-12*np.sin(2*np.pi*50*times)+1e-12*np.cos(2*np.pi*100*times)
    raw=mne.io.RawArray(values,mne.create_info(['MEG0111','MISC001','STI101'],sfreq,['mag','misc','stim']),
        first_samp=2500,verbose='ERROR')
    if date:raw.set_meas_date(datetime(2026,1,1,tzinfo=timezone.utc))
    raw.set_annotations(mne.Annotations([.4,1.+length/sfreq],[.6,.5-length/sfreq],['BAD_blink','BAD_muscle']))
    if length==0:raw.annotations.append(raw.first_time+.8,0.,'EDGE boundary')
    before=raw.copy();reference=raw.copy()
    # An independent public-Raw call provides a reference on every usable
    # span. The temporary mask is used only in this disposable test reference.
    if 0<length<=8:
        reference.annotations.append(reference.first_time+1.,length/sfreq,'BAD_reference_fragment')
    if length==0:
        reference.annotations.delete(np.flatnonzero(reference.annotations.description=='EDGE boundary'))
    reference.notch_filter(picks='mag',skip_by_annotation=('bad','edge'),verbose='ERROR',**SETTINGS)
    skipped=guarded_filter(raw)
    np.testing.assert_array_equal(raw.get_data(),reference.get_data())
    assert raw.first_samp==before.first_samp and raw.n_times==before.n_times and raw.info['sfreq']==before.info['sfreq']
    np.testing.assert_array_equal(raw.annotations.onset,before.annotations.onset)
    np.testing.assert_array_equal(raw.annotations.duration,before.annotations.duration)
    np.testing.assert_array_equal(raw.annotations.description,before.annotations.description)
    np.testing.assert_array_equal(raw.get_data()[1:],before.get_data()[1:])
    for start,stop in skipped:
        np.testing.assert_array_equal(raw.get_data(start=start,stop=stop),before.get_data(start=start,stop=stop))
    events=np.array([[raw.first_samp+1000,0,1],[raw.first_samp+5000,0,1]])
    kwargs=dict(event_id={'condition':1},tmin=-.5,tmax=2.667,baseline=None,detrend=1,
        proj=False,picks='all',reject={'mag':4e-12},reject_by_annotation=True,preload=True,verbose='ERROR')
    epochs=mne.Epochs(raw,events,**kwargs)
    # Restore annotations on the reference before epoch rejection, so the
    # comparison uses precisely the same accepted scientific mask.
    reference._annotations = before.annotations.copy()
    np.testing.assert_array_equal(reference.annotations.onset, raw.annotations.onset)
    other=mne.Epochs(reference,events,**kwargs)
    assert epochs.selection.tolist()==other.selection.tolist()==[1]
    np.testing.assert_array_equal(epochs.get_data(),other.get_data())


@pytest.mark.parametrize('length',list(range(1,33))+[63,64,127,128,255,256,999,1000,4999,5000,5001,9999,10000,10001,14999,15000,15001,20000,20001,30001])
def test_supported_lengths_against_real_cluster_notch(length):
    """Check both small DPSS limits and the long-window overlap/tail cases."""
    values=np.random.default_rng(8).normal(size=(1,length))*1e-13
    if length<=8:
        with pytest.raises((ValueError,ZeroDivisionError,AssertionError)):
            mne.filter.notch_filter(values,1000.,verbose='ERROR',**SETTINGS)
    else:
        out=mne.filter.notch_filter(values,1000.,verbose='ERROR',**SETTINGS)
        assert out.shape==values.shape and np.isfinite(out).all()
