"""Check that local read acceleration preserves bytes and source validation.

Usage: python -B -m pytest -q meg/tests/test_verified_read_cache.py
All files are tiny temporary fixtures; no participant sources are accessed.
"""
from pathlib import Path
import sys

import mne
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code/preprocessing_basic'))
from verified_read_cache import stage_read_cache, verified_fif_reads


# Exercise real FIF reads, then independent failure paths at the source/cache.
def test_verified_cache_preserves_samples_annotations_and_original_stamp(tmp_path):
    raw = mne.io.RawArray(np.arange(1000, dtype=float)[None] * 1e-12,
                         mne.create_info(['MEG0111'], 200., 'mag'), verbose='ERROR')
    raw.set_annotations(mne.Annotations([.4], [.1], ['BAD_blink']))
    source = tmp_path / 'source_raw.fif'
    raw.save(source, verbose='ERROR')
    reader = mne.io.read_raw_fif
    original = reader(source, preload=True, verbose='ERROR')
    certificate = stage_read_cache([source], tmp_path / 'cache')
    before = source.stat()
    with verified_fif_reads(certificate):
        cached = mne.io.read_raw_fif(source, preload=True, verbose='ERROR')
        np.testing.assert_array_equal(cached.get_data(), original.get_data())
        assert cached.annotations == original.annotations
        assert cached.first_samp == original.first_samp
    assert mne.io.read_raw_fif is reader
    assert source.stat().st_mtime_ns == before.st_mtime_ns


def test_verified_cache_refuses_changed_original_and_modified_cache(tmp_path):
    source = tmp_path / 'source_raw.fif'
    source.write_bytes(b'source bytes')
    certificate = stage_read_cache([source], tmp_path / 'cache')
    local = Path(certificate['sources'][0]['cache'])
    local.chmod(0o644)
    local.write_bytes(b'changed copy')
    local.chmod(0o444)
    with pytest.raises(RuntimeError, match='content changed'):
        with verified_fif_reads(certificate):
            pass
    other = stage_read_cache([source], tmp_path / 'cache2')
    source.write_bytes(b'new source bytes')
    with pytest.raises(RuntimeError, match='stamp or local content'):
        with verified_fif_reads(other):
            pass


def test_verified_cache_rechecks_source_during_reads_and_restores_reader(tmp_path):
    source = tmp_path / 'source_raw.fif'
    source.write_bytes(b'original')
    certificate = stage_read_cache([source], tmp_path / 'cache')
    reader = mne.io.read_raw_fif
    with pytest.raises(RuntimeError, match='changed during'):
        with verified_fif_reads(certificate):
            source.write_bytes(b'changed original')
            with pytest.raises(RuntimeError, match='before cached read'):
                mne.io.read_raw_fif(source)
    assert mne.io.read_raw_fif is reader
