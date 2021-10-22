# Sam(oa)²

Sam(oa)² stands for *Space-Filling Curves and Adaptive Meshes for Oceanic And Other Applications* and has been developed at TU Munich. Sam(oa)² is a PDE framework that can e.g., simulate wave propagations for tsunamis using adaptive meshes and space filling curves. For further information please refer to the original [README](README_samoa.md).

This repository contains the source code and scripts to build and run the following versions:

* A **work-sharing** version where in the main time stepping loop every thread in a process has a separate range of sections (work packages) assigned to process. This version might exhibit load imbalances between threads. Although that work distribution is comparable to a static schedule in an OpenMP work-sharing construct, the structure of the program does now easily allow to apply a dynamic schedule without severe code modifications.
* A **tasking** version that uses an over-decomposition approach with OpenMP tasks for the different work packages and allows a better work distribution between threads as idle threads might steal tasks from OpenMP task queues to mitigate the load imbalance.
* A **tasking** version that uses **Chameleon** to additionally apply task-based load balancing also between MPI ranks. Note that tasks in Chameleon are similar but not equal to OpenMP tasks.

## 0. Requirements

* Sam(oa)² requires the scons build system (see Prerequisites in the original README)
* In order to build the Chameleon version it is required to build Chameleon (https://github.com/chameleon-hpc/chameleon)
* As these experiments have been executed on the CLAIX supercomputer at RWTH Aachen University, the scripts might contain cluster specific commands/flags. Please feel free to adapt the scripts to your needs. Currently, the default environment uses an Intel compiler and Intel MPI

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

# compile the Chameleon version
# Append linux env vars with Chameleon include and lib folder (here: realized using an environment module)
module load chameleon
export CHAMELEON_DIR=/path/to/Chameleon/install/dir
COMPILE_CHAMELEON=1 ./samoa_build_intel.sh
```

## 2. Run

* In order to run the versions execute the following

```bash
# change to the scripts directoy
cd scripts/claix

# specify the paths to the Sam(oa)² and ASAGI main directory
export OUTPUT_DIR=/path/to/temp/output/folder # not used but required by program
export SAMOA_DIR=/path/to/samoa
export ASAGI_DATA_DIR=/path/to/asagi/dir

# run the work-sharing version
OMP_NUM_THREADS_VAR=11 NUM_RANKS=4 RUN_WS=1 ./samoa_run_intel.sh

# run the tasking version
OMP_NUM_THREADS_VAR=11 NUM_RANKS=4 RUN_TASKING=1 ./samoa_run_intel.sh

# run the Chameleon version
# Append linux env vars with Chameleon include and lib folder (here: realized using an environment module)
module load chameleon
OMP_NUM_THREADS_VAR=11 NUM_RANKS=4 RUN_CHAMELEON=1 ./samoa_run_intel.sh

# NOTE: it is also possible to execute the runs using a batch system like SLURM.
#       this will create a run with 4 processes per node
OMP_NUM_THREADS_VAR=11 NUM_RANKS=4 RUN_WS=1 RUN_TASKING=1 RUN_CHAMELEON=1 \
  sbatch --nodes=1 --ntasks-per-node=4 --cpus-per-task=12 \
  --export=OUTPUT_DIR,SAMOA_DIR,ASAGI_DATA_DIR,RUN_WS,RUN_TASKING,RUN_CHAMELEON,OMP_NUM_THREADS_VAR,NUM_RANKS \
  samoa_run_intel.sh
```

* Results of the run will be piped to log files

## 3. Traces

* The scripts additionally support building the versions with tracing tools like Extrae or Intel Trace Analyzer. For further information please refer to the script parameters