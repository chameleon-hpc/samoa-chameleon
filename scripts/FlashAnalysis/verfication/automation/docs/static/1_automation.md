# Automating the compilation, execution and rendering of the scenarios

The directory `scripts/FlashAnalysis/verification/automation` contains several bash-scripts to help with batch processing of the required scenarios.

All further directories will be specified relative to this path, unless specified otherwise.

Nearly all scripts log their output to `/logs`.

Each run is always executed for all limiters defined in `config.sh`.
 
## Requirements

### System packages

 - `bash`, with `tee`, `sed`, `grep`, `cat` and standard file handling tools like `mkdir` and `cp`
 - `pandoc` (the latest), `inkscape` and some kind of `pdflatex` or `xelatex` for the compilation to PDF
 - `python3`, `pip3`, `node`

### Node.js packages

 - `markdown-folder-to-html` npm module for the compilation to HTML

### Pip packages

 - `numpy` for computation
 - `scipy` for advanced interpolation
 - `numba` for JIT LLVM compilation
 - `pandas` for dataset processing
 - `psutil` for controlling worker threads
 - `matplotlib` for plotting
 - `pandocfilters` for the compilation to PDF

## Environment variables

 - `SCONSTHREADS`: Amount of parallel jobs for scons compilations
 - `THREADS`: Amount of OpenMP threads for execution
 - `RANKS`: Amount of MPI ranks for execution
 - `XDMFDIR`: Location of HDF5 and FoX libraries for compilation
 - `INTELRC`: Location of bash script to source to load intel compilers
 - `TIME`: Job execution time for SLURM

## Compilation

The directory `machines/` contains compilation and execution instructions for different architectures and compilers. Currently available are:
 - `docker-gnu`: GNU compiler in docker container (see `scripts/docker`)
 - `docker-intel`: Intel compiler in docker container (see `scripts/docker`)
 - `gnu`: GNU compiler on current system
 - `intel`: Intel compiler on current system

The file `config.sh` currently defines several scenarios to build binaries for.

The `compile-all.sh <machine-name>` and `compile-single <machine-name> <scenario-name>` scripts are used to batch invoke these configurations.

## Execution

Similar to compilation, the `run-all.sh <machine-name>` and `run-single.sh <machine-name> <run-name>` invokes the execution of the samoa binaries, as defined in the run configurations in the `runs/` directory.

### SLURM enqueue

The `enqueue-all.sh <machine-name>` and `enqueue-single.sh <machine-name> <run-name>` scripts can be used to invoke the SLURM batch scheduler instead of running the binary directly.

Do not forget to rename or copy `slurm.template.example` to `slurm.template` and adjust your email address. Also, you may have to provide absolute paths for your home directory instead of `~`. Make sure to create the directory where job logs will be placed.

## Rendering

Similar to the other scripts, the `render-all.sh` and `render-single.sh <render-name>` invokes the rendering / post-processing using the python scripts, as defined in the configurations contained in the `renders/` directory.

## Documentation

This documentation can be built using the `document-all.sh <mode>` script. It compiles the markdown files contained in `docs/` to HTML, and combines them with the images from `../results`. The completed documentation can the be found at `../documentation`. The same can be done for PDF.

The documentation compilation script generates the compilation and execution commands directly from the configuration files mentioned above.

## Examples

`./compile-all.sh machines/gnu`

`./compile-single.sh machines/gnu resting_lake`

`SCONSTHREADS=28 XDMFDIR=~/samoa_xdmf_libs/parallel INTELRC=/dev/null ./compile-all.sh machines/intel`

`./run-all.sh machines/docker-gnu`

`./run-single.sh machines/gnu runs/linear_beach-light.sh`

`RANKS=4 THREADS=28 TIME='00:15:00' ./enqueue-single.sh machines/intel runs/resting_lake.sh`

`./render-single.sh runs/linear_beach-light.sh`

`./document-all.sh pdf`