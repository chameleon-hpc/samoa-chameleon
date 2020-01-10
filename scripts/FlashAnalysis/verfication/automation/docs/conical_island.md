# Scenario 4.6: Flow around a conical island

## Notes

Experimental solutions for specific points in time are provided by https://nctr.pmel.noaa.gov/benchmark/Laboratory/Laboratory_ConicalIsland/index.html and converted to CSV here.

The refinement level `20` was chosen, which results in 2097152 cells, the same as in [1]: see section 4.6, paragraph three.

The dry/wet tolerance of `10e-3` and the simulation time of `20` were directly taken from [1].

The courant number of `0.3` was found using trial and error. The average computation time step size is about `tba` seconds, while the fixed time step size in [1] is `0.0025` seconds.

A XDMF output filter is used to limit the output to the region around the island. This is to reduce disk space usage. A second dataset contains the whole domain at every second.

For the incoming wave, a boundary condition is set with a time-dependent water height function.

## Compilation and Execution

## Evaluation

The water height at the gauge positions are sampled for each output step and plotted. These results can then be compared with [1], see figure 21 and 22.

Heightmap plots at two different times for each scenario are generated for the full domain. These results can then be compared with [1], see figure 23.
