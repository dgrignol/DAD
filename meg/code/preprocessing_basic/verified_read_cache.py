#!/usr/bin/env python3
"""Stage byte-verified, read-only FIF copies to accelerate slow shared-drive reads.

Usage from an ordered participant workflow (never changes accepted sources)::

    certificate = stage_read_cache(source_paths, Path('/tmp/NEW-participant-cache'))
    with verified_fif_reads(certificate):
        # Existing MNE readers/config/fingerprints still receive ORIGINAL paths.
        render_validated_participant_reports()

Each original is streamed once into a new local directory. SHA-256 is computed
during copying and independently checked on the local file; source size/mtime
must remain unchanged. Local copies are made read-only. Before use all cache
hashes and original stamps are checked. Every redirected open and context exit
rechecks the original stamp. Unlisted readers are unchanged. This is an I/O
optimization, not a new data source or a source/model-validation exception.

Outputs: local FIF copies and certificate.json with exact original paths,
stamps, byte hashes and local locations. Never write a source, modify raw sample
values, change fit preparation or bypass the caller's acceptance/fingerprint
checks. Stage copies for only one bounded participant operation at a time;
independent participants belong in separate processes/scheduler tasks.
"""
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import os

import mne

from continuous_cleaning import file_hash, file_stamp
from preprocessing_stage_utils import utc_now


# -----------------------------------------------------------------------------
# Copy once; retain the original identity and independently verify local bytes
# -----------------------------------------------------------------------------

def stage_read_cache(sources, directory):
    """Return a JSON certificate for NEW verified local copies of source FIFs."""
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    rows = []
    for source in map(lambda p: Path(p).resolve(), sources):
        target = directory / source.name
        if target.exists():
            raise RuntimeError('Cache source basenames must be unique.')
        before = file_stamp(source)
        print(f'Read cache: {source.name}; {before["size"]} bytes', flush=True)
        checksum = hashlib.sha256()
        with source.open('rb') as reader, target.open('xb') as writer:
            for block in iter(lambda: reader.read(8 * 1024 * 1024), b''):
                writer.write(block)
                checksum.update(block)
            writer.flush()
            os.fsync(writer.fileno())
        if file_stamp(source) != before or target.stat().st_size != before['size']:
            raise RuntimeError('Original source changed during staging.')
        expected = checksum.hexdigest()
        if file_hash(target) != expected:
            raise RuntimeError('Local copy differs from the bytes read from the original.')
        target.chmod(0o444)
        rows.append({'original': str(source), 'stamp': before,
                     'sha256': expected, 'cache': str(target)})
    certificate = {'generated_utc': utc_now(), 'sources': rows,
                   'scope': 'read-only byte-identical I/O cache; original source bindings retained'}
    with (directory / 'certificate.json').open('x') as handle:
        json.dump(certificate, handle, indent=2)
        handle.write('\n')
    return certificate


# -----------------------------------------------------------------------------
# Redirect ONLY exact certified opens, preserving all original validation calls
# -----------------------------------------------------------------------------

@contextmanager
def verified_fif_reads(certificate):
    """Temporarily route certified MNE FIF opens to byte-identical local copies.

    The caller must run this in a dedicated participant process without other
    threads opening FIFs. Readers are restored on success or any exception.
    """
    sources = {str(Path(row['original']).resolve()): row for row in certificate['sources']}
    if len(sources) != len(certificate['sources']):
        raise RuntimeError('Duplicate original source in read-cache certificate.')
    for row in sources.values():
        if file_stamp(row['original']) != row['stamp'] or file_hash(row['cache']) != row['sha256']:
            raise RuntimeError('Read-cache original stamp or local content changed.')
        if Path(row['cache']).stat().st_mode & 0o222:
            raise RuntimeError('Read-cache files must be read-only.')
    reader = mne.io.read_raw_fif

    def certified_reader(fname, *args, **kwargs):
        row = sources.get(str(Path(fname).resolve()))
        if row is None:
            return reader(fname, *args, **kwargs)
        if file_stamp(row['original']) != row['stamp']:
            raise RuntimeError('Original source changed before cached read.')
        return reader(row['cache'], *args, **kwargs)

    mne.io.read_raw_fif = certified_reader
    try:
        yield
    finally:
        mne.io.read_raw_fif = reader
        for row in sources.values():
            if file_stamp(row['original']) != row['stamp']:
                raise RuntimeError('Original source changed during cached reads.')
