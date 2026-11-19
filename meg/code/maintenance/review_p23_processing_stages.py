#!/usr/bin/env python3
"""Compare P23 rejected trials in acquisition, SSS and epoch preparation.

Usage::

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
      meg/code/maintenance/review_p23_processing_stages.py --workers 4 \
      --output meg/sandbox/p23_processing_stage_review_20261117

Read original acquisition FIFs (before offline preprocessing), accepted native
SSS/annotated FIFs and the completed noica-20261117-v2 trial decisions. Select
one representative rejected trial per acquisition by median number of failing
sensors, plus the weakest/strongest run-1 rejection and one retained control.
No trial is selected by the appearance of a plot or a desired causal result.

For each trial reconstruct the exact 50/100 Hz spectrum-fit operation over
complete BAD/EDGE-delimited stretches touching the display, then apply the
existing full-epoch linear detrend. Verify all reproduced amplitude-drop sensor
names against the saved decisions. Four-panel rows show the original continuous
signal, SSS continuous signal, line-cleaned continuous signal, and the epoch.
Separate epoch-comparison figures apply the SAME linear detrend to each stage.
Per-panel vertical scales are explicit; raw environmental fields can be much
larger than the SSS residual. Constant display centering never changes PTP.

Output: PNGs, paired numerical slices, JSON measurements and an HTML index in
a new review directory. No source, mask, threshold, scientific output or job is
changed. Independent acquisitions use separate processes, bounded by --workers;
numerical threads must be limited by the caller. Publish the index after join.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import html
import json
import multiprocessing
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
from mne.annotations import _annotations_starts_stops
import numpy as np
from scipy.signal import detrend

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"preprocessing_basic"))
from continuous_cleaning import annotation_rows,digest,file_stamp


# -----------------------------------------------------------------------------
# Read immutable sources and reproduce only the selected native windows.
# -----------------------------------------------------------------------------

def clean_context(raw, picks, left, stop):
    """Reproduce spectrum-fit using complete native good stretches, not crops."""
    data=raw.get_data(picks=picks,start=left,stop=stop)
    starts,stops=_annotations_starts_stops(raw,("bad","edge"),invert=True)
    for a,b in zip(starts,stops):
        if a>=stop or b<=left or b-a<9:continue
        segment=raw.get_data(picks=picks,start=int(a),stop=int(b))
        mne.filter.notch_filter(segment,raw.info["sfreq"],[50.,100.],method="spectrum_fit",
            filter_length="10s",copy=False,verbose="ERROR")
        lo,hi=max(left,int(a)),min(stop,int(b))
        data[:,lo-left:hi-left]=segment[:,lo-int(a):hi-int(a)]
    return data


def review_acquisition(task):
    """Produce paired raw-to-epoch evidence for one acquisition process."""
    source,trials,report,derivatives,output=task
    D=Path(derivatives);output=Path(output)
    original_path=D.parent/"data/sub23_090926/meg"/source
    annotated_path=D/"preprocessing/sub-23/ses-01/sss"/(Path(source).stem+"_ann.fif")
    print(f"P23 {source}: reading acquired and SSS headers",flush=True)
    original=mne.io.read_raw_fif(original_path,preload=False,verbose="ERROR")
    sss=mne.io.read_raw_fif(annotated_path,preload=False,verbose="ERROR")
    expected=next(r for r in report["preparation_contract"]["definition"]["preparation"]["inputs"] if r["source_filename"]==source)
    if digest(annotation_rows(sss))!=expected["annotations_sha256"] or [file_stamp(p) for p in sss.filenames]!=expected["parts"]:
        raise RuntimeError("Reviewed source or annotations changed")
    if (original.first_samp,original.n_times,original.info["sfreq"])!=(sss.first_samp,sss.n_times,sss.info["sfreq"]):
        raise RuntimeError("Acquired and SSS native sample frames differ")
    names=[sss.ch_names[i] for i in mne.pick_types(sss.info,meg=True,exclude=[])]
    kinds={n:sss.get_channel_types()[sss.ch_names.index(n)] for n in names}
    good=[n not in sss.info["bads"] for n in names]
    sfreq=sss.info["sfreq"];results=[]
    for trial in trials:
        anchor=int(float(trial["epoch_event_sample_raw"]))
        left=max(0,anchor-sss.first_samp-round(2.5*sfreq));stop=min(sss.n_times,anchor-sss.first_samp+round(4.5*sfreq)+1)
        a=anchor-sss.first_samp-round(.5*sfreq)-left;b=anchor-sss.first_samp+round(2.667*sfreq)+1-left
        acquired=original.get_data(picks=names,start=left,stop=stop)
        spatial=sss.get_data(picks=names,start=left,stop=stop)
        prepared=clean_context(sss,names,left,stop)
        epochs=[detrend(x[:,a:b],axis=-1,type="linear") for x in (acquired,spatial,prepared)]
        ptps=[np.ptp(x,axis=-1) for x in epochs]
        failed=[n for i,n in enumerate(names) if good[i] and ptps[-1][i]>report["amplitude_reject"][kinds[n]]]
        recorded=[n for n in trial["mne_drop_reasons"].split(",") if n.startswith("MEG")]
        if set(failed)!=set(recorded):raise RuntimeError(f"Drop reasons not reproduced: {trial['trial_key']}")
        shown=list(dict.fromkeys(["MEG0211","MEG1321",*recorded,"MEG1311","MEG0241"]))[:4]
        indices=[names.index(n) for n in shown]
        times=(np.arange(left,stop)+sss.first_samp-anchor)/sfreq
        epoch_times=np.arange(-round(.5*sfreq),round(2.667*sfreq)+1)/sfreq
        annotations=[(ann["onset"]-anchor/sfreq,ann["duration"],str(ann["description"])) for ann in sss.annotations
            if ann["onset"]<(sss.first_samp+stop)/sfreq and ann["onset"]+ann["duration"]>(sss.first_samp+left)/sfreq
            and str(ann["description"]).lower().startswith("bad")]
        title=f"P23 acquisition {trial['acquisition_run']} | behavioral run {trial['behavioral_run']}, trial {trial['behavioral_trial']} | {trial['condition_name']} | identity {int(float(trial['trial_sequence_identity']))}"
        columns=["Original acquisition: before offline preprocessing","After SSS: before line removal","After 50/100 Hz line removal","Epoch after linear detrend"]
        fig,axes=plt.subplots(4,4,figsize=(20,11),sharex="col")
        for row,index in enumerate(indices):
            for col,values in enumerate([acquired,spatial,prepared,epochs[-1]]):
                ax=axes[row,col]
                y=values[index]*1e12
                if col<3:
                    y=y-np.median(y[a:b]);ax.plot(times,y,color=["#66737b","#26769a","#2d6654"][col],lw=.65)
                    ax.axvspan(-.5,2.667,color="#7db7d8",alpha=.16)
                else:
                    y=y-(y.max()+y.min())*.5;ax.axhspan(-2,2,color="#b7d7c2",alpha=.2)
                    ax.plot(epoch_times,y,color="#b55226" if ptps[-1][index]>4e-12 else "#26384a",lw=.65)
                    ax.axhline(-2,color="#687a6e",ls="--",lw=.6);ax.axhline(2,color="#687a6e",ls="--",lw=.6)
                for start,duration,label in annotations:ax.axvspan(start,start+duration,color="#d96a75",alpha=.2)
                ax.axvline(0,color="#426c8e",lw=.7);ax.grid(alpha=.1)
                ax.set_xlim((-.5,2.667) if col==3 else (-2.5,4.5))
                if row==0:ax.set_title(columns[col],fontsize=10)
                if col==0:ax.set_ylabel(shown[row]+"\npT")
                if row==3:ax.set_xlabel("Seconds from corrected stimulus onset")
                stage=min(col,2)
                ax.text(.02,.95,f"Epoch PTP after detrend: {ptps[stage][index]*1e12:.2f} pT",transform=ax.transAxes,va="top",fontsize=8,
                    bbox=dict(facecolor="white",alpha=.7,edgecolor="none"))
        fig.suptitle(title+f"\nRecorded: {'retained control' if not recorded else str(len(recorded))+' failing sensors'} | each panel has its own labelled vertical scale",fontsize=12)
        fig.tight_layout(rect=(0,0,1,.92));name=f"run-{int(trial['acquisition_run']):02d}_sample-{anchor}"
        overview=output/(name+"_continuous_to_epoch.png");fig.savefig(overview,dpi=140);plt.close(fig)
        # Apply the same epoch operator to the acquired and SSS signals so that
        # changes due to SSS are distinguishable from epoch detrending itself.
        fig,axes=plt.subplots(4,3,figsize=(16,10),sharex=True)
        for row,index in enumerate(indices):
            for col,stage in enumerate(epochs):
                ax=axes[row,col];y=stage[index]*1e12;y=y-(y.max()+y.min())*.5
                ax.plot(epoch_times,y,color=["#66737b","#26769a","#b55226"][col],lw=.7)
                ax.axhspan(-2,2,color="#b7d7c2",alpha=.18)
                ax.axhline(-2,color="#687a6e",ls="--",lw=.6);ax.axhline(2,color="#687a6e",ls="--",lw=.6)
                ax.axvline(0,color="#426c8e",lw=.7);ax.grid(alpha=.1);ax.set_xlim(-.5,2.667)
                ax.set_title(f"{shown[row]} | PTP {ptps[col][index]*1e12:.2f} pT",fontsize=10,loc="left")
                if col==0:ax.set_ylabel("pT")
                if row==3:ax.set_xlabel("Seconds from corrected stimulus onset")
        fig.suptitle(title+"\nSame epoch and linear detrend at each stage: original acquisition | SSS | SSS + line removal",fontsize=11)
        fig.tight_layout(rect=(0,0,1,.92));paired=output/(name+"_matched_epoch_operator.png");fig.savefig(paired,dpi=140);plt.close(fig)
        np.savez_compressed(output/(name+"_paired_slices.npz"),channel_names=np.array(shown),times=times,
            acquired=acquired[indices],SSS=spatial[indices],line_cleaned=prepared[indices],epoch_times=epoch_times,
            acquired_epoch=epochs[0][indices],SSS_epoch=epochs[1][indices],prepared_epoch=epochs[2][indices])
        results.append({"trial_key":trial["trial_key"],"acquisition_run":int(trial["acquisition_run"]),"behavioral_run":trial["behavioral_run"],
            "behavioral_trial":trial["behavioral_trial"],"condition":trial["condition_name"],"identity":int(float(trial["trial_sequence_identity"])),
            "recorded_amplitude_channels":recorded,"reproduced_amplitude_channels":failed,"shown_channels":shown,
            "PTP_pT":{label:{n:float(v[names.index(n)]*1e12) for n in shown} for label,v in zip(["acquired_same_epoch_detrend","SSS_same_epoch_detrend","prepared_epoch"],ptps)},
            "source_acquired_parts":[file_stamp(p) for p in original.filenames],"source_annotated_parts":expected["parts"],
            "continuous_figure":overview.name,"epoch_figure":paired.name})
        print(f"P23 {source}, anchor {anchor}: reproduced {len(failed)} failing sensors",flush=True)
    original.close();sss.close();return results


# -----------------------------------------------------------------------------
# Deterministic trial selection, bounded participant reads, joined publication.
# -----------------------------------------------------------------------------

def main():
    """Generate paired evidence without applying the new no-amplitude policy."""
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--workers",type=int,default=4);p.add_argument("--output",type=Path,required=True)
    p.add_argument("--derivatives",type=Path,default=Path("/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives"));args=p.parse_args()
    root=args.derivatives/"processing_runs/noica-20261117-v2"
    report=json.loads((root/"reports/sub-23/ses-01/04_epoching_report.json").read_text())
    with (root/"analysis/sub-23/ses-01/meg/sub-23_ses-01_task-Attention_run-all_trial_decisions.csv").open() as f:rows=list(csv.DictReader(f))
    selected=[]
    for run in range(1,6):
        failed=sorted([r for r in rows if int(r["acquisition_run"])==run and "MEG" in r["mne_drop_reasons"]],key=lambda r:len(r["mne_drop_reasons"].split(",")))
        selected.append(failed[len(failed)//2])
        if run==1:
            selected.extend([failed[0],failed[-1],next(r for r in rows if r["acquisition_run"]=="1" and r["included_final"]=="True")])
    args.output.mkdir(parents=True,exist_ok=False);results=[]
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context("spawn")) as pool:
        futures=[pool.submit(review_acquisition,(name,[r for r in selected if r["source_fif"]==name],report,str(args.derivatives),str(args.output))) for name in sorted({r["source_fif"] for r in selected})]
        for future in as_completed(futures):results.extend(future.result())
    results.sort(key=lambda r:(r["acquisition_run"],r["trial_key"]))
    (args.output/"measurements.json").write_text(json.dumps({"purpose":"source_vs_processing_review_only","trials":results},indent=2))
    cards="".join(f'<h2>Acquisition {r["acquisition_run"]}, behavioral run {r["behavioral_run"]}, trial {r["behavioral_trial"]}</h2><p>{html.escape(r["trial_key"])}: {len(r["recorded_amplitude_channels"])} recorded amplitude failures.</p><img src="{r["continuous_figure"]}"><img src="{r["epoch_figure"]}">' for r in results)
    (args.output/"index.html").write_text('<!doctype html><meta charset="utf-8"><title>P23 source-to-epoch review</title><style>body{font:17px system-ui;max-width:1750px;margin:30px auto;padding:0 20px;line-height:1.5}img{width:100%}p{max-width:1100px}</style><h1>P23: original acquisition through epoch preparation</h1><p>Read-only review of eight deterministically selected trials spanning every acquisition. The original acquisition is before offline SSS. Blue shading marks the full epoch; red marks existing BAD annotations. Panel scales differ and are labelled. The paired epoch plots use the same linear detrend at every processing stage. The green 4 pT band is the final post-SSS epoch threshold, shown at earlier stages for reference only: raw acquisition amplitudes are not the production rejection criterion. SSS changes spatial mixtures, so individual sensor waveforms need not remain identical. No trial rejection, annotation or source recording was changed for these figures.</p>'+cards)
    print("P23 review saved:",args.output/"index.html",flush=True)


if __name__=="__main__":main()
