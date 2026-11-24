"""Validate the noise metric, unit handling and analysis output isolation.

Usage: python -m pytest meg/tests/test_drsa_noise_normalization.py -q
Synthetic correlated MEG-like data exercise unit conversions and rank loss;
no recordings or production derivatives are needed or modified.
"""
from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'code/preprocessing_basic'))
from drsa_noise_normalization import fit_transform, gaussian_smooth
from drsa_analysis_variant import configure_variant, source_config, assert_result_normalization


def noise():
    """Generate mixed-unit, rank-deficient correlated data across four runs."""
    rng=np.random.default_rng(41)
    latent=rng.normal(size=(40,5,25));mix=rng.normal(size=(8,5))
    data=np.einsum('cf,eft->ect',mix,latent)*np.array([1e-13]*3+[1e-11]*5)[None,:,None]
    return data,np.repeat(np.arange(4),10)


def transformed(data, fitted):
    return np.einsum('dc,ect->edt',fitted['white'],data/fitted['sigma'][None,:,None])


def test_unit_conversion_and_sensor_permutation_preserve_metric():
    data,runs=noise();fit,_=fit_transform(data,runs)
    gains=np.array([1e15]*3+[1e13]*5)
    converted,_=fit_transform(data*gains[None,:,None],runs)
    np.testing.assert_allclose(transformed(data,fit),transformed(data*gains[None,:,None],converted),atol=2e-12)
    order=np.array([7,1,5,0,3,6,2,4]);permuted,_=fit_transform(data[:,order],runs)
    np.testing.assert_allclose(transformed(data,fit)[:,order],transformed(data[:,order],permuted),atol=2e-12)
    assert fit['rank']==5


def test_covariance_whitening_retains_only_supported_noise_directions():
    data,runs=noise();fit,_=fit_transform(data,runs)
    y=transformed(data,fit).transpose(0,2,1).reshape(-1,8);y-=y.mean(0)
    eigenvalues=np.linalg.eigvalsh(y.T@y/len(y))
    np.testing.assert_allclose(eigenvalues[:3],0,atol=1e-12)
    np.testing.assert_allclose(eigenvalues[3:],1,atol=.01)
    np.testing.assert_allclose(gaussian_smooth(np.ones((2,3,60)),15),1,atol=1e-15)


def test_default_namespace_and_explicit_none_do_not_redirect_inputs(monkeypatch,tmp_path):
    monkeypatch.delenv('DAD_DRSA_NORMALIZATION',raising=False)
    monkeypatch.delenv('DAD_DRSA_ANALYSIS_LABEL',raising=False)
    root=tmp_path/'processing_runs'/'v2'
    cfg={'layout':{'subject_bids_id':'sub-04','session_bids_id':'ses-01'},'matlab_drsa':{},
         'resolved_paths':{'processing_run_root':str(root),'derivatives_analysis_dir':str(root/'analysis/sub-04/ses-01'),
                           'derivatives_reports_dir':str(root/'reports/sub-04/ses-01')}}
    mnn=configure_variant(cfg);none=configure_variant(cfg,'none')
    assert mnn['matlab_drsa']['neural_normalization']=='mnn'
    assert mnn['resolved_paths']['derivatives_analysis_dir']!=none['resolved_paths']['derivatives_analysis_dir']
    assert source_config(mnn)['resolved_paths']==cfg['resolved_paths']
    assert not Path(mnn['resolved_paths']['derivatives_analysis_dir']).is_relative_to(root)
    with pytest.raises(ValueError,match='Mixed dRSA normalization'):
        assert_result_normalization({'neural_normalization':'none'},mnn)
