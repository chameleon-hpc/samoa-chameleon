# Sam(oa)²

* TODO: Intro
* TODO: Link to original Readme File
* Requirement: Sam(oa)² uses the scons build system
* As these experiments have been executed on the CLAIX supercomputer at RWTH Aachen University, the scripts might contain cluster specific commands. Currently, the default environment uses an Intel compiler and Intel MPI

This repository contains the source code and scripts to build and run the following 2 versions:

* A **work-sharing** version where in the main time stepping loop every thread in a process has a separate range of sections (work packages) assigned to process. This version might exhibit load imbalances between threads. Although that work distribution is comparable to a static schedule in an OpenMP work-sharing construct, the structure of the program does now easily allow to apply a dynamic schedule without severe code modifications.
* A **tasking** version that uses an over-decomposition approach with OpenMP tasks for the different work packages and allows a better work distribution between threads as idle threads might steal tasks from OpenMP task queues to mitigate the load imbalance.

## 1. Build

* In order to build the versions execute the following:

```bash
# change to the scripts directoy
cd scripts/claix

# specify the paths to the Sam(oa)² and ASAGI main directory
export SAMOA_DIR=/path/to/samoa
export ASAGI_DIR=/path/to/ASAGI

# compile the work-sharing version
COMPILE_WS=1 ./samoa_build_intel.sh

# compile the tasking version
COMPILE_TASKING=1 ./samoa_build_intel.sh
```

## 2. Run

* In order to run the versions execute the following

```bash
# change to the scripts directoy
cd scripts/claix

# specify the paths to the Sam(oa)² and ASAGI main directory
export SAMOA_DIR=/path/to/samoa
export OUTPUT_DIR=/path/to/temp/output

# run the work-sharing version
OMP_NUM_THREADS=11 NUM_RANKS=2 RUN_WS=1 ./samoa_run_intel.sh

# run the tasking version
OMP_NUM_THREADS=11 NUM_RANKS=2 RUN_TASKING=1 ./samoa_run_intel.sh

# NOTE: it is also possible to execute the runs using a batch system like SLURM
OMP_NUM_THREADS=11 NUM_RANKS=2 RUN_WS=1 RUN_TASKING=1 \
  sbatch --export=OUTPUT_DIR,SAMOA_DIR,RUN_WS,RUN_TASKING,OMP_NUM_THREADS,NUM_RANKS \
  samoa_run_intel.sh
```

* Results of the run will be piped to log files