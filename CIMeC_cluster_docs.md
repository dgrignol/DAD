### *Legacy notes: CIMeC SGE cluster*

> **Status:** This guide describes the previous CIMeC cluster and its Sun Grid
> Engine (`qsub`) scheduler. For current access, connect to
> `hpc.dipsen.unitn.it`; the DAD project directory remains
> `/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD` for now. The old
> host, direct compute-node SSH, and `qsub` instructions below are historical
> and must not be used on DIPSEN. See the [current DIPSEN HPC documentation](https://hpc.dipsen.unitn.it/)
> for login and SLURM instructions.

The remaining sections are historical examples for accessing the old cluster,
submitting parallel jobs with `qsub`, and running code on individual nodes.
A commented submission script is provided at the end for reference only.

---

#### *1\. Access the Cluster*

##### **1.1 Connect to the VPN (if outside UNITN network)**

Follow instructions on the IT ServiceDesk page

(UNITN VPN: [https://servicedesk.unitn.it/sd/it/service/vpn?id=unitrento\_v2\_service\_card\&table=sc\_cat\_item\&sys\_id=e5a898fec35bfd104cbb7055df013171\&searchTerm=vpn](https://servicedesk.unitn.it/sd/it/service/vpn?id=unitrento_v2_service_card&table=sc_cat_item&sys_id=e5a898fec35bfd104cbb7055df013171&searchTerm=vpn) ).

##### **1.2 Open a Terminal**

* macOS: **Terminal**  
* Windows: **Windows Terminal** (you can do it from PowerShell too but commands may differ slightly)

Run:

```
ssh -Y <your.unitn.username>@cluster.cimec.unitn.it
```

Example:

```
(base) damiano@MacBook ~ % ssh -Y damiano.grignolio@cluster.cimec.unitn.it
```

Enter your password (not visible while typing).

After login you will see something like:

```
[damiano.grignolio@cluster ~]$
```

You are now in the **HEAD NODE**.

**Do not run any script or computation here.**

This node is only for file management and dispatching jobs.

---

#### *2\. Running Code Directly on a Compute Node (useful for debugging)*

You can log into a specific compute node (e.g. compute-0-3):

```
ssh compute-0-3.local

```

(You can check node names and load here:

http://cluster.cimec.unitn.it/ganglia/?r=hour\&m=load\_one\&c=cluster )

Once on the node, enter your project directory

example

```
cd /mnt/storage/tier2/morwur/Projects/DAMIANO
```

##### **2.1 Activate your virtual environment**

If your code requires python3 or specific libraries, you **must** activate the correct virtual environment e.g. a micromamba environment (if you do not know what I am talking about see BONUS 7\. Create Micromamba Environment): 

example

```
micromamba activate drsa3111
(drsa311) [damiano.grignolio@compute-0-3 DAMIANO]$
```

You can set up a micromamba environment from the terminal of your computer (I advise you to set it up in your home folder on the cluster so you can always activate it from anywhere else in the cluster. Ask chat GPT how to set up one with the packages you require). 

##### **2.2 Run a script manually**

Example:

```
python dRSA_run.py --subjects 1-14
```

This will run sequentially on that node only.

---

#### *3\. Submitting Jobs with qsub (run in parallel)*

Navigate to your project folder: Run “ssh \-Y \<your.unitn.username\>@[cluster.cimec.unitn.it](http://cluster.cimec.unitn.it)” from the terminal and navigate to your project directory. All projects live under /mnt/storage/tier2/morwur/ ....

Example:

```
cd /mnt/storage/tier2/morwur/Projects/DAMIANO
```

Your prompt becomes:

```
[damiano.grignolio@cluster DAMIANO]$
```

You submit the bash (.sh) script from the **head node**:

```
[name.surname@cluster DAMIANO]$ qsub <your_script.sh>
```

The .sh file defines:

* the job name  
* resources requested  
* environment activation  
* which Python/MATLAB script to run  
* arguments (subject ID, model, session, etc.)

See below on how to specify the bash script.

---

#### *4\. Minimal Job Submission Script (qsub)*

Create an .sh file. Example:

```
submit_example.sh
```

##### **Example: submit\_example.sh** 

```
#!/bin/bash
#$ -N dRSA_analysis                  # Give a name to the job
#$ -S /bin/bash                # Use bash
#$ -cwd                        # Run in current working dir
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.out   # Stdout log per task
#$ -e logs/$JOB_NAME.$JOB_ID.$TASK_ID.err   # Stderr log per task
#$ -j y                        # Merge stdout/stderr (useful has you will end up with just one file for logs)
#$ -pe smp 4                   # Request 4 CPU cores
#$ -l h_vmem=24G               # Request 24 GB RAM
#$ -t 1-27                     # Array: one task per subject (1..27)


set -euo pipefail              # Exit on errors and undefined vars


WD="/mnt/storage/tier2/morwur/Projects/DAMIANO/"  # Project path
mkdir -p logs                  # Ensure log folder exists
cd "$WD"                       # Work inside the project


# Activate the environment (micromamba example)
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$($HOME/bin/micromamba shell hook -s bash)" # or wherever it sits
micromamba activate drsa311
export PYTHONNOUSERSITE=1      # Avoid user-site packages


# Pick subject from the array index (01, 02, …)
subject=$(printf "%02d" "$SGE_TASK_ID")
subject_label="sub-${subject}"


echo ">>> Running my_script_dRSA for ${subject_label}"
python my_script_dRSA.py "${subject_label}"
```

Key knobs you may change:

\- \`-t 1-27\` to set how many array tasks (e.g., number of subjects).

\- \`-pe smp 4\` for CPU cores; match to your script.

\- \`-l h\_vmem=24G\` for memory per task.

\- The Python command at the bottom (add flags as needed).

For running matlab code you don’t need an environment (like micromamba above) but you need to call matlab (without UI). Example:

```
#!/bin/bash

### SGE Options ###
# 
# Job name 
#$ -N BioMotion
# Script interpreter
#$ -S /bin/bash  
# Log file
#$ -o my_matlab_job.log
# Merge error and output
#$ -j yes 

###################

# Define arguments to submitted Matlab function
WD=$1
XYZ=$2
iSub=$3

# Call Matlab
/state/partition1/MATLAB/R2020a/bin/matlab -nodesktop -nosplash -nojvm -nodisplay -r "s3_startSVM_on_cluster('$WD','$XYZ',$iSub)"

exit
```

#### *5\. Submit the Example Job*

Submit with optional script arguments (after \`--\` to separate from \`qsub\` options): 

examples

```
qsub submit_example.sh
qsub submit_example.sh -- --t 1-3 # Run first 3 sub only
# for the matlab example above qsub -wd $WD ./s2_submit_matlab_function.sh $WD xyz 16
```

---

#### *6\. Monitoring Jobs*

Run:

```
qstat
```

\- \`qstat\` shows running/pending jobs with their job-ID and node.

\- \`qstat \-j \<job-ID\>\` gives details; \`qdel \<job-ID\>\` cancels.

\- Logs live under \`logs/\` with the names set in the script.

Example output:

```
job-ID  prior  name     user            state  start_at                queue                             slots
87563   0.50   C1_dRSA  damiano.grig    r      11/19/2025 18:03:07     all.q@compute-0-1.local          4
88138   0.50   C1_dRSA  damiano.grig    r      11/26/2025 11:42:52     all.q@compute-0-3.local          4
```

All console output appears in the log files generated by qsub:

```
<jobname>.o<jobID>   # STDOUT
<jobname>.e<jobID>   # STDERR
```

---

#### *BONUS 7\. Create Micromamba Environment*

Minimal guide: Python3 environment on Cluster using micromamba

##### **7.1 Install micromamba (once)**

Navigate to your home directory as we will install the environment there.

```shell
cd ~
mkdir -p bin
cd bin
curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba --strip-components=1
./micromamba shell init -s bash -r $HOME/micromamba
exec $SHELL
```

##### **7.2 Create the environment**

Here we name the environment “drsa311” but you can change it to your likings.

```shell
micromamba create -y -n drsa311 python=3.11 pip
micromamba activate drsa311
```

##### **7.3 Install packages** 

(I give some example packages but you can choose yours)

```shell
micromamba install -y -c conda-forge \
  numpy pandas scipy scikit-learn \
  matplotlib seaborn \
  mne mne-bids nibabel \
  pysoundfile openpyxl tqdm wordfreq
```

##### **7.4 Test the environment**

```shell
micromamba activate drsa311
python -c "import numpy, mne; print('OK')"
```

You can use the environment from any project folder

The environment is installed in your home directory (\~/micromamba) and can be used from any location.

Example

```shell
cd /path/to/your/project
micromamba activate drsa311
python pipeline_wrapper.py \
  --subjects 2-23 
```

Use the environment inside cluster jobs (SGE)

Each job runs in a fresh shell, so you must activate the environment inside the job script.

```shell
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$($HOME/bin/micromamba shell hook -s bash)"
micromamba activate drsa311
```

##### **Key points**

The environment lives in: \~/micromamba/envs/drsa311

It works from any folder.

It must be activated in every new shell and inside every cluster job script.

You can create as many independent environments as you want just by creating another with a different name.

---

#### *Quick checklist*

\- VPN on (if remote) → SSH to head.

\- Choose: \`cd\` project → interactive node (\`ssh compute-0-x.local\`) **OR** \`qsub\` from head.

\- Activate the env environment before running or in the bash script.

\- For \`qsub\`, set job name, resources, array range, and log paths; check \`qstat\`.

Each subject becomes an independent job → the cluster scheduler spreads them across nodes.

