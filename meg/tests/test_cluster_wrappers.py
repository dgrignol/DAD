"""
test_cluster_wrappers.py

Purpose
Validate the SLURM submission wrappers without submitting jobs or running
scientific preprocessing. Verify that the apply wrapper delegates through
DAD_PROJECT_ROOT and that P24 Stage 1 stops after the digital trigger stage
with the deep report enabled; timing adapters belong to a later cohort gate.

Usage
Run only these tests:
  python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_cluster_wrappers.py -q
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import shlex
import subprocess
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
APPLY_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_preprocessing_apply.sh"
DRSA_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_drsa_analysis.sh"
DRSA_GROUP_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_drsa_group_average.sh"
DRSA_LINE_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_drsa_line_summaries.sh"
DOT_SIDE_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_dot_side_decoding.sh"
DOT_SIDE_GROUP_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_dot_side_group_report.sh"
DOT_SIDE_TG_WRAPPER = (
    REPO_ROOT / "meg" / "cluster" / "submit_dot_side_time_generalization.sh"
)
DOT_SIDE_TG_GROUP_WRAPPER = (
    REPO_ROOT
    / "meg"
    / "cluster"
    / "submit_dot_side_time_generalization_group_report.sh"
)
DEPLOY_WRAPPER = REPO_ROOT / "meg" / "cluster" / "deploy_review_to_cluster.sh"
STAGE1_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_preprocessing_stage1.sh"
STAGE2_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_preprocessing_stage2.sh"


@pytest.mark.parametrize("mode", ["render", "recount", "failed_recount", "missing_input", "bad_label", "array", "existing_visual"])
def test_retention_html_wrapper_renders_saved_counts_or_recounts_first(tmp_path, mode):
    """Dispatch only requested counting/rendering stages and stop on failure."""
    project=tmp_path/"project"; project.mkdir()
    derivatives=tmp_path/"derivatives"
    reports=derivatives/"annotation_triage/reports"
    baseline=reports/"before"; baseline.mkdir(parents=True)
    (baseline/"report.json").write_text("fixture")
    if mode not in ("recount","failed_recount","missing_input"):
        current=reports/"now"; current.mkdir()
        (current/"report.json").write_text("fixture")
    binary=tmp_path/"bin"; binary.mkdir()
    mamba=binary/"micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n");mamba.chmod(0o755)
    python=binary/"python"
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$DAD_TEST_OUTPUT"\n'
        'if [[ "${DAD_TEST_FAIL_COUNT:-0}" == "1" && "$*" == *annotation_triage.py* ]]; then exit 9; fi\n')
    python.chmod(0o755)
    captured=tmp_path/"calls.txt"
    env=dict(os.environ,PATH=str(binary)+os.pathsep+os.environ["PATH"],
        DAD_PROJECT_ROOT=str(project),DAD_DERIVATIVES_ROOT=str(derivatives),
        MICROMAMBA_BIN=str(mamba),DAD_TRIAGE_LABEL="now",DAD_TRIAGE_VISUAL_LABEL="html-now",
        DAD_TRIAGE_COMPARE_TO="before",DAD_TRIAGE_RECOUNT="1" if mode in ("recount","failed_recount") else "0",
        DAD_TEST_OUTPUT=str(captured),SLURM_CPUS_PER_TASK="4",TMPDIR=str(tmp_path/"jobtmp"))
    env.pop("SLURM_ARRAY_TASK_ID",None)
    if mode=="failed_recount":env["DAD_TEST_FAIL_COUNT"]="1"
    if mode=="bad_label":env["DAD_TRIAGE_LABEL"]="../other"
    if mode=="array":env["SLURM_ARRAY_TASK_ID"]="4"
    if mode=="existing_visual":(derivatives/"annotation_triage/visual_reports/html-now").mkdir(parents=True)
    wrapper=REPO_ROOT/"meg/cluster/submit_annotation_retention_report.sh"
    result=subprocess.run(["bash",str(wrapper)],env=env,capture_output=True,text=True)
    calls=captured.read_text().splitlines() if captured.exists() else []
    if mode not in ("render","recount"):
        assert result.returncode!=0
        assert len(calls)==(1 if mode=="failed_recount" else 0)
        assert not any("render_annotation_triage_report.py" in c for c in calls)
        return
    assert result.returncode==0,result.stdout+result.stderr
    assert len(calls)==(2 if mode=="recount" else 1)
    assert shlex.split(calls[-1])==["-B",str(project/"meg/code/preprocessing_basic/render_annotation_triage_report.py"),
        "--report-dir",str(reports/"now"),"--label","html-now","--jobs","4","--compare-to",str(reports/"before")]
    if mode=="recount":
        tokens=shlex.split(calls[0]); assert "--subjects" in tokens
        assert tokens[tokens.index("--subjects")+1:tokens.index("--label")]==list(map(str,range(4,40)))
    assert not any("import_ocular" in c or "run_pipeline.py" in c for c in calls)
    for name in ("submit_annotation_retention_report.sh","render_annotation_triage_report.py"):
        assert name in DEPLOY_WRAPPER.read_text()


@pytest.mark.parametrize("failure", [None, "import", "array", "existing_snapshot", "missing_notes"])
def test_annotation_triage_validates_import_before_cohort_counting(tmp_path, failure):
    """Verify ordered cohort dispatch and fail-closed behavior without data/jobs."""
    project = tmp_path / "project"; project.mkdir()
    derivatives = tmp_path / "derivatives"
    (derivatives / "annotation_triage/component_reports/test-books").mkdir(parents=True)
    decisions = tmp_path / "completed.csv"; decisions.write_text("fixture\n")
    binary = tmp_path / "bin"; binary.mkdir()
    mamba = binary / "micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n")
    mamba.chmod(0o755)
    runtime = binary / "python"
    runtime.write_text(
        '#!/bin/bash\nprintf "%s\\n" "$*" >> "$DAD_TEST_OUTPUT"\n'
        'if [[ "${DAD_TEST_FAIL_IMPORT:-0}" == "1" && "$*" == *import_ocular_component_choices.py* ]]; then exit 9; fi\n'
    )
    runtime.chmod(0o755)
    captured = tmp_path / "calls.txt"
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
        DAD_PROJECT_ROOT=str(project), DAD_DERIVATIVES_ROOT=str(derivatives),
        MICROMAMBA_BIN=str(mamba), DAD_COMPONENT_REPORT_LABEL="test-books",
        DAD_TRIAGE_LABEL="test-retention", DAD_TRIAGE_DECISIONS=str(decisions),
        DAD_TEST_OUTPUT=str(captured), SLURM_CPUS_PER_TASK="4", TMPDIR=str(tmp_path / "jobtmp"))
    env.pop("SLURM_ARRAY_TASK_ID", None)
    if failure == "import":
        env["DAD_TEST_FAIL_IMPORT"] = "1"
    elif failure == "array":
        env["SLURM_ARRAY_TASK_ID"] = "15"
    elif failure == "existing_snapshot":
        (derivatives / "annotation_triage/reports/test-retention").mkdir(parents=True)
    elif failure == "missing_notes":
        decisions.unlink()
    wrapper = REPO_ROOT / "meg/cluster/submit_annotation_triage.sh"
    result = subprocess.run(["bash", str(wrapper)], env=env, capture_output=True, text=True)
    calls = captured.read_text().splitlines() if captured.exists() else []
    if failure:
        assert result.returncode != 0
        assert len(calls) == (1 if failure == "import" else 0)
        assert not any("annotation_triage.py" in call for call in calls)
        return
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(calls) == 2
    config_args = ["--config", str(project / "meg/params/pipeline_default.yaml"),
        "--site-config", str(project / "meg/params/pipeline_cluster.yaml"),
        "--subject-config-dir", str(project / "meg/params")]
    assert shlex.split(calls[0]) == [
        "-B", str(project / "meg/code/preprocessing_basic/import_ocular_component_choices.py"),
        *config_args, "--reports-dir", str(derivatives / "annotation_triage/component_reports/test-books"),
        "--decisions", str(decisions), "--reviewer", "Damiano", "--apply", "--jobs", "4"]
    assert shlex.split(calls[1]) == [
        "-B", str(project / "meg/code/preprocessing_basic/annotation_triage.py"),
        *config_args, "--subjects", *map(str, range(4, 40)), "--label", "test-retention", "--jobs", "4"]
    assert not any("run_pipeline.py" in call or "sbatch" in call for call in calls)
    assert "meg/cluster/submit_annotation_triage.sh" in DEPLOY_WRAPPER.read_text()


def test_component_report_wrapper_selects_only_offline_book_generation(tmp_path):
    """Dispatch through harmless runtimes; no scheduler/scientific command may run."""
    project = tmp_path / "project"; project.mkdir()
    binary = tmp_path / "bin"; binary.mkdir()
    mamba = binary / "micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n")
    mamba.chmod(0o755)
    runtime = binary / "python"
    runtime.write_text('#!/bin/bash\nprintf "%s\\n" "$*" > "$DAD_TEST_OUTPUT"\n')
    runtime.chmod(0o755)
    captured = tmp_path / "calls.txt"
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
        DAD_PROJECT_ROOT=str(project), MICROMAMBA_BIN=str(mamba), DAD_COMPONENT_REPORT_LABEL="test-books",
        SLURM_ARRAY_TASK_ID="15", DAD_TEST_OUTPUT=str(captured), TMPDIR=str(tmp_path / "jobtmp"))
    wrapper = REPO_ROOT / "meg/cluster/submit_ocular_component_reports.sh"
    result = subprocess.run(["bash", str(wrapper)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    command = captured.read_text()
    assert shlex.split(command) == [
        "-B", str(project / "meg/code/preprocessing_basic/build_ocular_component_reports.py"),
        "--config", str(project / "meg/params/pipeline_default.yaml"),
        "--site-config", str(project / "meg/params/pipeline_cluster.yaml"),
        "--subject-config-dir", str(project / "meg/params"),
        "--subjects", "15", "--label", "test-books", "--resume",
    ]
    assert "run_pipeline.py" not in command and "sbatch" not in command
    deployed = DEPLOY_WRAPPER.read_text()
    for name in ("build_ocular_component_reports.py", "import_ocular_component_choices.py",
                 "ocular_component_report.html", "ocular_component_report.js", "submit_ocular_component_reports.sh"):
        assert name in deployed


@pytest.mark.parametrize("stage,rebuild", [("artifact_annotation", "0"),
                                         ("artifact_annotation", "1"), ("ica_manual", "1")])
def test_affected_cleaning_refresh_commands_select_only_the_requested_stage(tmp_path, stage, rebuild):
    """Test user command construction with harmless runtimes; submit no jobs."""
    project = tmp_path / "project"
    params = project / "meg/params"; params.mkdir(parents=True)
    for name in ("pipeline_default.yaml", "pipeline_cluster.yaml", "sub-04.yaml"):
        (params / name).write_text("fixture: true\n")
    calibration = project / "meg/meg_calibration_files/fwmegmaxfilter"; calibration.mkdir(parents=True)
    for name in ("ct_sparse_triux2.fif", "sss_cal_3182_250410.dat"):
        (calibration / name).write_bytes(b"fixture")
    report = tmp_path / "derivatives/reports/sub-04/ses-01"; report.mkdir(parents=True)
    (report / "03_tuning_choices.yaml").write_text("fixture: true\n")
    binary = tmp_path / "bin"; binary.mkdir()
    mamba = binary / "micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n"); mamba.chmod(0o755)
    python = binary / "python"
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$DAD_TEST_OUTPUT"\n'); python.chmod(0o755)
    captured = tmp_path / "calls.txt"
    processing_pass = "apply" if stage == "ica_manual" else "review"
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
        DAD_PROJECT_ROOT=str(project), DAD_DERIVATIVES_ROOT=str(tmp_path / "derivatives"),
        MICROMAMBA_BIN=str(mamba), MAMBA_ROOT_PREFIX=str(tmp_path / "mamba"),
        TMPDIR=str(tmp_path / "jobtmp"), SLURM_ARRAY_TASK_ID="4", DAD_TEST_OUTPUT=str(captured),
        DAD_MAXWELL_QC_ONLY="0", DAD_BLINK_AID_ONLY="0", DAD_STAGES=stage,
        DAD_PROCESSING_PASS=processing_pass, DAD_REBUILD=rebuild)
    result = subprocess.run(["bash", str(REPO_ROOT / "meg/cluster/submit_preprocessing_review.sh")],
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    calls = captured.read_text().splitlines()
    assert len(calls) == 2  # Runtime imports, then exactly one selected runner.
    assert f"--stages {stage} --processing-pass {processing_pass}" in calls[1]
    assert ("--overwrite" if rebuild == "1" else "--resume") in calls[1]
    assert "--tuning-config" in calls[1]
    assert "maxwell_filter" not in calls[1] and "trigger_extraction" not in calls[1]


def test_p24_stage1_runs_only_digital_trigger_stage(tmp_path):
    """Run the actual wrapper against harmless environment/Python stubs.

    P24 must make exactly one runner call and leave deep reports enabled.
    No scheduler, actual pipeline, raw input, or derivative tree is accessed.
    A reintroduced P24 adapter invocation would produce a second stub call.
    """
    project = tmp_path / "project"
    params = project / "meg" / "params"
    params.mkdir(parents=True)
    for name in ("pipeline_default.yaml", "pipeline_cluster.yaml", "sub-24.yaml"):
        (params / name).write_text("fixture: true\n", encoding="utf-8")
    binary_dir = tmp_path / "bin"
    binary_dir.mkdir()
    mamba = binary_dir / "micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n", encoding="utf-8")
    mamba.chmod(0o755)
    python = binary_dir / "python"
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$DAD_TEST_OUTPUT"\n', encoding="utf-8")
    python.chmod(0o755)
    captured = tmp_path / "runner_calls.txt"
    env = dict(os.environ, PATH=str(binary_dir) + os.pathsep + os.environ["PATH"],
               DAD_PROJECT_ROOT=str(project), MICROMAMBA_BIN=str(mamba),
               MAMBA_ROOT_PREFIX=str(tmp_path / "mamba"), TMPDIR=str(tmp_path / "jobtmp"),
               SLURM_ARRAY_TASK_ID="24", DAD_TEST_OUTPUT=str(captured))
    result = subprocess.run(["bash", str(STAGE1_WRAPPER)], env=env,
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    calls = captured.read_text(encoding="utf-8").splitlines()
    assert len(calls) == 1
    assert "run_pipeline.py" in calls[0]
    assert "--stages trigger_extraction --write-inventory" in calls[0]
    assert "--disable-deep-trigger-report" not in calls[0]
    assert "p24_acquisition_adapter.py" not in calls[0]


def test_maxwell_plot_refresh_never_invokes_preprocessing_runner(tmp_path):
    """Exercise the actual plot-only wrapper path using harmless runtime stubs."""
    project = tmp_path / "project"
    params = project / "meg/params";params.mkdir(parents=True)
    for name in ["pipeline_default.yaml", "pipeline_cluster.yaml", "sub-04.yaml"]:
        (params / name).write_text("fixture: true\n")
    calibration = project / "meg/meg_calibration_files/fwmegmaxfilter";calibration.mkdir(parents=True)
    for name in ["ct_sparse_triux2.fif", "sss_cal_3182_250410.dat"]:
        (calibration / name).write_bytes(b"fixture")
    binary = tmp_path / "bin";binary.mkdir()
    mamba = binary / "micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n");mamba.chmod(0o755)
    python = binary / "python"
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$DAD_TEST_OUTPUT"\n');python.chmod(0o755)
    captured = tmp_path / "calls.txt"
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
        DAD_PROJECT_ROOT=str(project), DAD_DERIVATIVES_ROOT=str(tmp_path / "derivatives"),
        DAD_RAW_ROOT=str(tmp_path / "raw"), MICROMAMBA_BIN=str(mamba),
        MAMBA_ROOT_PREFIX=str(tmp_path / "mamba"), TMPDIR=str(tmp_path / "jobtmp"),
        SLURM_ARRAY_TASK_ID="4", DAD_TEST_OUTPUT=str(captured), DAD_MAXWELL_QC_ONLY="1",
        DAD_QC_WINDOW_SEC="20", DAD_QC_SEED="97", DAD_PROCESSING_PASS="review",
        DAD_REBUILD="0", DAD_STAGES="maxwell_filter")
    wrapper = REPO_ROOT / "meg/cluster/submit_preprocessing_review.sh"
    result = subprocess.run(["bash", str(wrapper)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    calls = captured.read_text().splitlines()
    assert len(calls) == 2  # Runtime imports, followed by the plot-only CLI.
    assert "maxwell_qc.py" in calls[1] and "--window-sec 20 --seed 97 --overwrite" in calls[1]
    assert "--subjects 04" in calls[1]
    assert not any("run_pipeline.py" in call for call in calls)


@pytest.mark.parametrize("rebuild", ["0", "1"])
@pytest.mark.parametrize("aid", ["blink", "squid"])
def test_preparation_wrapper_does_not_run_sss_muscle_or_final_ica(tmp_path, rebuild, aid):
    """Execute preparation command construction with dummy runtimes, no jobs."""
    project = tmp_path / "project"
    params = project / "meg/params"; params.mkdir(parents=True)
    for name in ("pipeline_default.yaml", "pipeline_cluster.yaml", "sub-04.yaml"):
        (params / name).write_text("fixture: true\n")
    calibration = project / "meg/meg_calibration_files/fwmegmaxfilter"; calibration.mkdir(parents=True)
    for name in ("ct_sparse_triux2.fif", "sss_cal_3182_250410.dat"):
        (calibration / name).write_bytes(b"fixture")
    report = tmp_path / "derivatives/reports/sub-04/ses-01"; report.mkdir(parents=True)
    (report / "03_tuning_choices.yaml").write_text("fixture: true\n")
    binary = tmp_path / "bin"; binary.mkdir()
    mamba = binary / "micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n"); mamba.chmod(0o755)
    python = binary / "python"
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$DAD_TEST_OUTPUT"\n'); python.chmod(0o755)
    captured = tmp_path / "calls.txt"
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
        DAD_PROJECT_ROOT=str(project), DAD_DERIVATIVES_ROOT=str(tmp_path / "derivatives"),
        MICROMAMBA_BIN=str(mamba), MAMBA_ROOT_PREFIX=str(tmp_path / "mamba"),
        TMPDIR=str(tmp_path / "jobtmp"), SLURM_ARRAY_TASK_ID="4", DAD_TEST_OUTPUT=str(captured),
        DAD_MAXWELL_QC_ONLY="0", DAD_BLINK_AID_ONLY=str(int(aid == "blink")),
        DAD_SQUID_PREPARE_ONLY=str(int(aid == "squid")), DAD_PROCESSING_PASS="review", DAD_REBUILD=rebuild,
        DAD_SOURCE="sub04_run01_220526.fif")
    env.pop("DAD_STAGES", None)  # Blink selects 4B; SQUID keeps the 4A default.
    result = subprocess.run(["bash", str(REPO_ROOT / "meg/cluster/submit_preprocessing_review.sh")],
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    calls = captured.read_text().splitlines()
    assert len(calls) == 2
    entry = "blink_annotation.py" if aid == "blink" else "review_squid_jumps.py"
    assert entry in calls[1] and "--source sub04_run01_220526.fif" in calls[1]
    assert ("--overwrite" in calls[1]) == (aid == "blink" and rebuild == "1")
    assert ("--resume" in calls[1]) == (aid == "blink" and rebuild == "0")
    if aid == "squid":
        assert "--mode prepare" in calls[1]
    assert all(name not in calls[1] for name in ("run_pipeline.py", "03_A_MaxFilt.py", "03_B_ArtefactAnnotation.py", "03_C_ICA.py"))
    assert entry in DEPLOY_WRAPPER.read_text()


@pytest.mark.parametrize("mode", ["subject", "summary"])
def test_stage2_wrapper_runs_only_trial_preparation(tmp_path, mode):
    """Real wrappers invoke one standalone accounting operation, with no runner."""
    project = tmp_path / "project"
    cluster = project / "meg/cluster"
    cluster.mkdir(parents=True)
    shutil.copy2(STAGE2_WRAPPER, cluster / STAGE2_WRAPPER.name)
    binary = tmp_path / "bin"
    binary.mkdir()
    mamba = binary / "micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n")
    mamba.chmod(0o755)
    python = binary / "python"
    python.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$DAD_TEST_OUTPUT"\n')
    python.chmod(0o755)
    captured = tmp_path / "calls"
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
               DAD_PROJECT_ROOT=str(project), MICROMAMBA_BIN=str(mamba),
               MAMBA_ROOT_PREFIX=str(tmp_path / "mamba"), TMPDIR=str(tmp_path / "jobtmp"),
               SLURM_ARRAY_TASK_ID="30", DAD_TEST_OUTPUT=str(captured))
    for key in ("DAD_STAGE2_MODE", "DAD_STAGE2_OVERWRITE", "DAD_STAGE2_VALIDATE_ONLY"):
        env.pop(key, None)
    wrapper = STAGE2_WRAPPER if mode == "subject" else STAGE2_WRAPPER.with_name("submit_preprocessing_stage2_summary.sh")
    result = subprocess.run(["bash", str(wrapper)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    calls = captured.read_text().splitlines()
    assert len(calls) == 1
    assert "trial_preparation.py" in calls[0]
    assert "run_pipeline.py" not in calls[0]
    assert "--overwrite" not in calls[0]
    assert ("--summarize-cohort" in calls[0]) == (mode == "summary")
    assert ("sub-30.yaml" in calls[0]) == (mode == "subject")


def test_stage2_deployment_includes_the_complete_accounting_contract():
    """The new standalone stage, both wrappers and technical guide ship together."""
    content = DEPLOY_WRAPPER.read_text()
    for name in ("trial_preparation.py", "submit_preprocessing_stage2.sh",
                 "submit_preprocessing_stage2_summary.sh", "PIPELINE_SPEC_basic.md"):
        assert name in content


def test_apply_wrapper_resolves_review_wrapper_from_project_root(tmp_path):
    """
    A copied apply script must delegate through DAD_PROJECT_ROOT.

    This test starts the apply wrapper from an unrelated directory and uses a
    harmless review-wrapper stub to record the delegated environment.
    """
    project_root = tmp_path / "project"
    cluster_dir = project_root / "meg" / "cluster"
    cluster_dir.mkdir(parents=True)
    output_path = tmp_path / "delegated.txt"

    review_wrapper = cluster_dir / "submit_preprocessing_review.sh"
    review_wrapper.write_text(
        "#!/bin/bash\n"
        'printf "%s|%s\\n" "$DAD_PROCESSING_PASS" "$DAD_REBUILD" '
        '>"$DAD_TEST_OUTPUT"\n',
        encoding="utf-8",
    )
    review_wrapper.chmod(0o755)

    external_dir = tmp_path / "external_copy"
    external_dir.mkdir(parents=True)
    copied_apply = external_dir / "submit_preprocessing_apply.sh"
    shutil.copy2(APPLY_WRAPPER, copied_apply)
    copied_apply.chmod(0o755)

    env = os.environ.copy()
    env.update(
        {
            "DAD_PROJECT_ROOT": str(project_root),
            "DAD_TEST_OUTPUT": str(output_path),
        }
    )
    completed = subprocess.run(
        [str(copied_apply)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert completed.returncode == 0, completed.stderr
    assert output_path.read_text(encoding="utf-8") == "apply|1\n"


def test_drsa_and_deploy_wrappers_have_valid_bash_syntax():
    """Cluster submission/deployment shell scripts must parse before upload."""
    for script in (
        DRSA_WRAPPER,
        DRSA_GROUP_WRAPPER,
        DRSA_LINE_WRAPPER,
        DOT_SIDE_WRAPPER,
        DOT_SIDE_GROUP_WRAPPER,
        DOT_SIDE_TG_WRAPPER,
        DOT_SIDE_TG_GROUP_WRAPPER,
        DEPLOY_WRAPPER,
        *sorted((REPO_ROOT / "meg" / "cluster").glob("submit_*.sh")),
    ):
        completed = subprocess.run(
            ["bash", "-n", str(script)],
            check=False,
            capture_output=True,
            text=True,
        )
        assert completed.returncode == 0, completed.stderr


def test_deploy_payload_includes_matlab_analysis_layer():
    """Stage-07 deployment must include MATLAB scripts and helper functions."""
    text = DEPLOY_WRAPPER.read_text(encoding="utf-8")
    assert '"meg/code/analysis/"' in text
    assert '"meg/code/preprocessing_basic/07_run_matlab_drsa.py"' in text
    assert '"meg/code/preprocessing_basic/08_plot_drsa_line_summaries.py"' in text
    assert '"meg/code/preprocessing_basic/drsa_result_utils.py"' in text
    assert '"meg/cluster/submit_drsa_analysis.sh"' in text
    assert '"meg/cluster/submit_drsa_group_average.sh"' in text
    assert '"meg/cluster/submit_drsa_line_summaries.sh"' in text
    assert '"simulations/functions/"' in text


def test_deploy_payload_includes_dot_side_sanity_check():
    """Optional dot-side cluster jobs need their script and wrapper deployed."""
    text = DEPLOY_WRAPPER.read_text(encoding="utf-8")
    assert '"meg/code/sanity_checks/"' in text
    assert '"meg/cluster/submit_dot_side_decoding.sh"' in text
    assert '"meg/cluster/submit_dot_side_group_report.sh"' in text
    assert '"meg/cluster/submit_dot_side_time_generalization.sh"' in text
    assert '"meg/cluster/submit_dot_side_time_generalization_group_report.sh"' in text


def test_deploy_payload_includes_dot_erf_wrapper():
    """The optional dot-ERF job wrapper must ship with its Python scripts."""
    text = DEPLOY_WRAPPER.read_text(encoding="utf-8")
    assert '"meg/cluster/submit_dot_erf_analysis.sh"' in text


def test_job_wrappers_use_slurm_directives_and_runtime_variables():
    """Every submit wrapper should use SLURM names, not SGE scheduler state."""
    for script in sorted((REPO_ROOT / "meg" / "cluster").glob("submit_*.sh")):
        text = script.read_text(encoding="utf-8")
        assert "#SBATCH --job-name=" in text, script.name
        assert "#SBATCH --output=" in text, script.name
        assert "#SBATCH --error=" in text, script.name
        assert "#SBATCH --time=" in text, script.name
        assert "#$ -" not in text, script.name
        assert "SGE_TASK_ID" not in text, script.name
        assert "NSLOTS" not in text, script.name
        assert "qsub" not in text, script.name


def test_dot_side_wrapper_exposes_decoding_controls():
    """The dot-side wrapper should expose the requested runtime knobs."""
    text = DOT_SIDE_WRAPPER.read_text(encoding="utf-8")
    assert "DAD_DOT_SIDE_ANALYSES" in text
    assert "DAD_DOT_SIDE_CONDITION" in text
    assert "DAD_DOT_SIDE_MODEL_SOURCE" in text
    assert "occluded_deviant_predicted" in text
    assert "DAD_DECODE_WINDOW_SAMPLES" in text
    assert "DAD_TRIALS_PER_AVERAGE" in text
    assert "DAD_RANDOM_SEED" in text
    assert "--validate-only" in text


def test_dot_side_group_wrapper_exposes_cluster_controls():
    """The group wrapper should expose CI/SEM and cluster-threshold settings."""
    text = DOT_SIDE_GROUP_WRAPPER.read_text(encoding="utf-8")
    assert "DAD_DOT_SIDE_ANALYSES" in text
    assert "--analyses" in text
    assert "DAD_DOT_SIDE_INTERVAL" in text
    assert "DAD_DOT_SIDE_SMOOTH_SIGMA_SAMPLES" in text
    assert "--smooth-sigma-samples" in text
    assert "DAD_CLUSTER_FORMING_P" in text
    assert "DAD_CLUSTER_ALPHA" in text
    assert "DAD_DOT_SIDE_GROUP_SUBJECTS" in text
    assert "submit_dot_side_group_report.sh" in text


def test_dot_side_time_generalization_wrappers_expose_cluster_controls():
    """Time-generalization wrappers should expose matrix and 2D cluster knobs."""
    subject_text = DOT_SIDE_TG_WRAPPER.read_text(encoding="utf-8")
    group_text = DOT_SIDE_TG_GROUP_WRAPPER.read_text(encoding="utf-8")

    assert "03_dot_side_time_generalization.py" in subject_text
    assert "DAD_DOT_SIDE_TG_ANALYSES" in subject_text
    assert "DAD_DECODE_WINDOW_SAMPLES" in subject_text
    assert "DAD_TRIALS_PER_AVERAGE" in subject_text
    assert "all_trials_per_identity" in subject_text
    assert "dot_side_time_generalization" in subject_text

    assert "04_dot_side_time_generalization_group_report.py" in group_text
    assert "DAD_DOT_SIDE_TG_METRIC" in group_text
    assert "DAD_DOT_SIDE_TG_SMOOTH_SIGMA_SAMPLES" in group_text
    assert "DAD_DOT_SIDE_SMOOTH_SIGMA_SAMPLES" in group_text
    assert "--smooth-sigma-samples" in group_text
    assert "DAD_DOT_SIDE_TG_INCLUDE_CLUSTER_TABLE" in group_text
    assert "--include-cluster-table" in group_text
    assert "DAD_CLUSTER_FORMING_P" in group_text
    assert "DAD_CLUSTER_ALPHA" in group_text
    assert "DAD_CLUSTER_TAIL" in group_text
    assert "2D cluster" in group_text


def test_drsa_line_summary_wrapper_exposes_stage08_controls():
    """The Stage-08 wrapper should expose group line-summary controls."""
    text = DRSA_LINE_WRAPPER.read_text(encoding="utf-8")
    assert "DAD_DRSA_LINE_SUBJECTS" in text
    assert "DAD_DRSA_LINE_LAG_WINDOW_SEC" in text
    assert "DAD_DRSA_LINE_INCLUDE_SMOOTHED" in text
    assert "DAD_DRSA_LINE_SHADING" in text
    assert "DAD_DRSA_LINE_INTERVAL_LAG_MODE" in text
    assert "DAD_DRSA_LINE_MODEL_TIME_WINDOW_SAMPLES" not in text
    assert "--model-time-window-samples" not in text
    assert 'SHADING="${DAD_DRSA_LINE_SHADING:-across-subject-sem}"' in text
    assert "--shading" in text
    assert "--interval-lag-mode" in text
    assert "08_plot_drsa_line_summaries.py" in text


def test_drsa_group_average_wrapper_exposes_stage07_controls():
    """The Stage-07 group-average wrapper should expose useful cluster controls."""
    text = DRSA_GROUP_WRAPPER.read_text(encoding="utf-8")
    assert "DAD_DRSA_GROUP_CONFIG" in text
    assert "DAD_DRSA_GROUP_SITE_CONFIG" in text
    assert "DAD_DRSA_GROUP_SUBJECT_CONFIG" in text
    assert "DAD_DRSA_GROUP_SUBJECTS" in text
    assert "DAD_DRSA_GROUP_SMOOTHING_PLOTS" in text
    assert "DAD_DRSA_GROUP_PLOT_SMOOTHING_WINDOW" in text
    assert "--group-average-only" in text
    assert "--group-subjects" in text
    assert "--print-summary" in text
    assert "07_run_matlab_drsa.py" in text


def test_drsa_wrapper_exposes_model_source_controls():
    """Stage-07 cluster submissions should expose observed/predicted model runs."""
    text = DRSA_WRAPPER.read_text(encoding="utf-8")
    assert "DAD_DRSA_MODEL_SOURCES" in text
    assert "DAD_DRSA_CONDITIONS" in text
    assert "pe_prior_constrained_50" in text


def test_drsa_wrapper_uses_dipsen_matlab_module_and_memory():
    """Stage 07 should discover the DIPSEN MATLAB module and retain JVM memory."""
    text = DRSA_WRAPPER.read_text(encoding="utf-8")
    assert "module load matlab/R2020a" in text
    assert "module load matlab/R2025b" in text
    assert "#SBATCH --mem=24G" in text
