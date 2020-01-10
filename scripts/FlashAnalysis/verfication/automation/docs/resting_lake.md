# Scenario 4.1: Lake at rest

## Notes

Two simulation data sets with different temporal resolutions are provided, one with `400` and one with `4000` output steps. For evaluation, the 4000-step variant is used, the 400-step variant is used only for development and testing purposes.

The refinement level `11` was chosen, which results in a triangle leg length of about `0.022` and `4096` cells - the same as in [1]: see section 4.1, paragraph two, and compare with figure 4, left panel.

The dry/wet tolerance of `10e-6` and simulation time of `40` seconds were directly taken from [1].

The courant number of `0.3` was found using trial and error. The average computation time step size is about `0.001333` resp. `0.00125` seconds, while the fixed time step size in [1] is `0.002` seconds.

## Compilation and Execution

## Evaluation

The L2 and Lsup limiter errors are computed and subsequently plotted by comparing the numerical results to the analytical solution. These results can then compared with [1], see figure 3 and 4.
