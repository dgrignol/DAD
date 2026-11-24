"""Resolve isolated sensor-noise analysis variants for Stages 07--10.

Usage: pass --neural-normalization mnn (default) or none to any dRSA stage.
DAD_DRSA_NORMALIZATION provides the same cluster-wide selection; optional
DAD_DRSA_ANALYSIS_LABEL chooses a fresh simple label. Inputs remain the
accepted Stage-04/06 processing run. Outputs, reports and manifests go to a
sibling processing-run directory, outside legacy recursive result searches.
This module never changes global preprocessing paths or existing derivatives.
"""
from copy import deepcopy
import os
from pathlib import Path
import re
from pipeline_config import load_and_resolve_config as _load


def add_analysis_arguments(parser):
    """Expose one explicit normalization selector consistently in all stages."""
    parser.add_argument('--neural-normalization', choices=('mnn', 'none'), default=None,
                        help='Sensor-noise normalization; default mnn. none reproduces the legacy metric.')


def source_config(cfg):
    """Recover immutable epoch/export paths when validating upstream lineage."""
    source = deepcopy(cfg)
    if 'drsa_input_resolved_paths' in source:
        source['resolved_paths'] = source.pop('drsa_input_resolved_paths')
    return source


def configure_variant(cfg, mode=None):
    """Bind one mode and output namespace without redirecting its input files."""
    if 'drsa_input_resolved_paths' in cfg:
        return cfg
    cfg = deepcopy(cfg)
    mode = mode or os.environ.get('DAD_DRSA_NORMALIZATION') or cfg.get('matlab_drsa', {}).get('neural_normalization', 'mnn')
    if mode not in ('mnn', 'none'):
        raise ValueError('neural_normalization must be mnn or none.')
    cfg.setdefault('matlab_drsa', {})['neural_normalization'] = mode
    resolved = cfg['resolved_paths']
    if mode == 'none' and not resolved.get('processing_run_root') and not os.environ.get('DAD_DRSA_ANALYSIS_LABEL'):
        return cfg  # Explicit access to historical unlabelled results.
    cfg['drsa_input_resolved_paths'] = deepcopy(resolved)
    source_root = Path(resolved.get('processing_run_root', Path(resolved['derivatives_analysis_dir']).parents[2]))
    label = os.environ.get('DAD_DRSA_ANALYSIS_LABEL', f'{source_root.name}-drsa-{mode}-v1')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,119}', label):
        raise ValueError('DAD_DRSA_ANALYSIS_LABEL must be a simple directory label.')
    target = (source_root.parent / label if resolved.get("processing_run_root") else source_root / "analysis_variants" / label)
    if target.resolve() == source_root.resolve():
        raise ValueError('dRSA output namespace must differ from its input processing run.')
    cfg['drsa_analysis_root'] = str(target)
    cfg['drsa_analysis_label'] = label
    for key, directory in [('derivatives_analysis_dir', 'analysis'), ('derivatives_reports_dir', 'reports')]:
        resolved[key] = str(target / directory / cfg['layout']['subject_bids_id'] / cfg['layout']['session_bids_id'])
    return cfg


def load_drsa_config(*args, neural_normalization=None, **kwargs):
    """Load normal pipeline configuration, then select analysis-only outputs."""
    cfg, warnings, config_path, subject_path = _load(*args, **kwargs)
    return configure_variant(cfg, neural_normalization), warnings, config_path, subject_path


def assert_result_normalization(metadata, cfg):
    """Refuse mixed neural metrics before downstream grouping or inference."""
    if 'drsa_input_resolved_paths' not in cfg:
        return
    expected = cfg['matlab_drsa']['neural_normalization']
    actual = metadata.get('neural_normalization', 'none')
    if actual != expected:
        raise ValueError(f'Mixed dRSA normalization: expected {expected}, found {actual}.')
