# Scenario 4.5: Runup onto a complex three-dimensional beach

## Notes

Experimental solutions for specific points in time are provided by http://isec.nacse.org/workshop/2004_cornell/bmark2.html and converted to CSV here.

The refinement level `18` was chosen, which results in tba cells, which is roughly the same for the rectangular section where the bathymetry is given as in [1]: see section 4.6, paragraph three.

The dry/wet tolerance of `10e-4` and the simulation time of `40` were directly taken from [1].

The courant number of `0.3` was found using trial and error. The average computation time step size is about `tba` seconds, while the fixed time step size in [1] is `0.001` seconds.

A XDMF output filter is used to limit the output to a rectangle at the shore which contains the gauges.

To fill the square domain, the bathymetry of the most bottom position was extended downwards.

For the incoming wave, a boundary condition is set with time-dependent water height values read from a file.

## Compilation and Execution

## Evaluation

The water height at the gauge positions are sampled for each output step and plotted. These results can then be compared with [1], see figure 18.

Heightmap plots at `t = 15.0, 15.5, 16.0, 16.5, 17.0` are generated for the full domain. These results can then be compared with [1], see figure 19.