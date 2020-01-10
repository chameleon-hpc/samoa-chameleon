# Scenario 4.2: Tsunami runup onto a linearly sloping beach

## Notes

Analytical solutions for specific points in time are provided by http://isec.nacse.org/workshop/2004_cornell/bmark1.html and converted to CSV here.

The refinement level `20` was chosen, which results in a triangle leg length of about `49.219‬` - roughly the same as in [1] (`50`): see section 4.2, paragraph two.

The dry/wet tolerance of `10e-2` was directly taken from [1]. The simulation time of `230` seconds was chosen to fit all interesting points in time.

The courant number of `0.35` was found using trial and error. The average computation time step size is about `0.0155` seconds, while the fixed time step size in [1] is `0.04` seconds.

A XDMF output filter is used to limit the output to a small horizontal slice in the middle. This is to reduce disk space usage. Furthermore, unused HDF5 chunks were removed manually from the dataset (The `-light` dataset).

The y axis offset was changed such that the bathymetry is zero at the coordinate origin. This was done to ensure the correct calculation of the water level offset.

## Compilation and Execution

## Evaluation

A z-cross-section is interpolated and plotted for all relevant variables, and compared to the analytical solution at `t = 160, 175, 220`. These results can then be compared with [1], see figure 6.

Line `19` in `linear_beach.py` may be modified to interpolate at different z values (default: `25000`).
