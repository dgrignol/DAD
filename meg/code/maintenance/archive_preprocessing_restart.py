#!/usr/bin/env python3
"""Archive the approved pilot and first cohort Stage-1 attempt in place.

Purpose
-------
Preserve pilot outputs for P01-P07 and their group analyses under
``derivatives/pilot``. Preserve the first P08-P39 Stage-1 attempt separately
under ``derivatives/previous_runs/stage1_2026-10-15``. This leaves the active
preprocessing/report locations ready for a fresh P04-P39 run.

Usage from the DAD repository root
--------------------------------
Preview the exact source/destination list (the default is read-only)::

    python meg/code/maintenance/archive_preprocessing_restart.py

Apply the user-approved archive once::

    python meg/code/maintenance/archive_preprocessing_restart.py --apply

An alternate derivative root can be supplied with ``--derivatives-root``.
Only existing derivative directories in the explicit plan are considered.
All sources are inventoried before the first move. Directory renames on the
same filesystem preserve file contents and internal relative paths. Each
destination is checked against its original path/size/mtime inventory, with
the inventory digest and byte/file counts recorded in an archive manifest.
The manifest's digest describes metadata, not a content checksum.

Safety and provenance
---------------------
The sibling ``Attention/data`` tree is never opened or modified. The script
rejects paths under a directory named ``data``, symlinked archive roots,
destination collisions, and files modified within the last five minutes
(Finder metadata is exempt). Existing JSON/YAML provenance remains verbatim;
the archive manifest records the mapping from historical to archived paths.
An interrupted move is recorded in the manifest and requires inspection;
rerunning this one-time command cannot overwrite the archive.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import time
from datetime import datetime, timezone


DEFAULT_ROOT = Path("/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives")


def inventory_tree(root: Path) -> list[dict]:
    """Record files, symlinks, and empty directories without reading contents."""
    entries = []
    for directory, subdirs, files in os.walk(root, followlinks=False):
        folder = Path(directory)
        if not subdirs and not files:
            entries.append({"path": folder.relative_to(root).as_posix(), "kind": "empty_directory"})
        for name in sorted(subdirs + files):
            path = folder / name
            info = path.lstat()
            if stat.S_ISLNK(info.st_mode):
                entries.append({"path": path.relative_to(root).as_posix(),
                                "kind": "symlink", "target": os.readlink(path)})
            elif stat.S_ISREG(info.st_mode):
                entries.append({"path": path.relative_to(root).as_posix(),
                                "kind": "file", "size_bytes": info.st_size,
                                "mtime_ns": info.st_mtime_ns})
    return sorted(entries, key=lambda item: (item["path"], item["kind"]))


def approved_plan(root: Path) -> list[tuple[Path, Path, Path]]:
    """Return exact approved moves with each associated archive root."""
    pilot = root / "pilot"
    previous = root / "previous_runs/stage1_2026-10-15"
    plan = []
    for category in ("preprocessing", "reports"):
        for number in range(1, 40):
            relative = Path(category) / f"sub-{number:02d}"
            archive = pilot if number <= 7 else previous
            plan.append((root / relative, archive / relative, archive))
    for relative in ("reports/group", "analysis", "sanity_checks", "poster_figures"):
        plan.append((root / relative, pilot / relative, pilot))
    return plan


def write_manifest(path: Path, payload: dict) -> None:
    """Atomically checkpoint the archive journal outside all moved sources."""
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> int:
    """Inventory the fixed plan, optionally rename, then verify every move."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--derivatives-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    root = args.derivatives_root.expanduser().resolve(strict=True)
    if "data" in root.parts or root.name != "derivatives":
        raise ValueError(f"Expected a derivatives root outside raw data: {root}")
    plan = [(source, target, archive) for source, target, archive in approved_plan(root)
            if source.exists()]
    for source, target, archive in plan:
        if source.is_symlink() or not source.is_dir():
            raise ValueError(f"Expected an ordinary source directory: {source}")
        if target.exists() or target.is_symlink():
            raise FileExistsError(f"Preserving existing archive destination: {target}")
        if archive.is_symlink() or not archive.resolve().is_relative_to(root):
            raise ValueError(f"Unsafe archive location: {archive}")
        print(f"{source.relative_to(root)} -> {target.relative_to(root)}", flush=True)
    if not args.apply:
        print(f"Preview only: {len(plan)} directory moves; pass --apply after approval.")
        return 0

    # Inventory every source first. Reject recent scientific output activity;
    # the user's completed jobs and explicit approval authorize this migration.
    snapshots = {}
    for source, _, _ in plan:
        entries = inventory_tree(source)
        if any(item["kind"] == "file" and Path(item["path"]).name != ".DS_Store"
               and item["mtime_ns"] / 1e9 > time.time() - 300 for item in entries):
            raise RuntimeError(f"Recent output activity; archive deferred: {source}")
        snapshots[source] = entries
        print(f"Inventoried {source.relative_to(root)}: {len(entries)} entries", flush=True)

    journals = {}
    for archive in sorted({archive for _, _, archive in plan}):
        manifest = archive / "archive_manifest.json"
        if manifest.exists():
            raise FileExistsError(f"Preserving existing archive journal: {manifest}")
        archive.mkdir(parents=True, exist_ok=True)
        journals[archive] = {"created_utc": datetime.now(timezone.utc).isoformat(),
                             "status": "moving", "verification": "relative path, type, size, mtime",
                             "historical_provenance": "Preserved verbatim; embedded paths may refer to original locations.",
                             "moves": []}
        with (archive / "archive_files.jsonl").open("x", encoding="utf-8") as handle:
            for source, target, assigned in plan:
                if assigned == archive:
                    for entry in snapshots[source]:
                        handle.write(json.dumps({"source_directory": str(source),
                                                 "destination_directory": str(target), **entry}) + "\n")
        write_manifest(manifest, journals[archive])

    # Renaming avoids copying large FIF/MAT contents. Verify each destination
    # against the full pre-move inventory before marking that move complete.
    for source, target, archive in plan:
        entries = snapshots[source]
        record = {"source": str(source), "destination": str(target), "status": "moving",
                  "file_count": sum(item["kind"] == "file" for item in entries),
                  "total_bytes": sum(item.get("size_bytes", 0) for item in entries),
                  "metadata_inventory_sha256": hashlib.sha256(
                      json.dumps(entries, sort_keys=True).encode()).hexdigest()}
        journal = journals[archive]
        journal["moves"].append(record)
        write_manifest(archive / "archive_manifest.json", journal)
        target.parent.mkdir(parents=True, exist_ok=True)
        source.rename(target)
        if inventory_tree(target) != entries:
            record["status"] = "verification_failed"
            write_manifest(archive / "archive_manifest.json", journal)
            raise RuntimeError(f"Archive metadata changed during move: {target}")
        record["status"] = "verified"
        write_manifest(archive / "archive_manifest.json", journal)
        print(f"Verified {target.relative_to(root)}: {record['file_count']} files, {record['total_bytes']} bytes", flush=True)
    for archive, journal in journals.items():
        journal["status"] = "complete"
        write_manifest(archive / "archive_manifest.json", journal)
    print(f"Archive complete: {len(plan)} verified directory moves.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
